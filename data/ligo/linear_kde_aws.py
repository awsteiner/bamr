#!/usr/bin/python3

import numpy as np
import getopt
import sys
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV
import h5py
import os

compute_bandwidth=False
sample=False
sample_file=''
plots=False
grid_file='gw170817_kde.o2'

try:
    opts, args=getopt.getopt(sys.argv[1:],'bghps:')
except getopt.GetoptError as err:
    print(str(err))
    quit()
for opt, arg in opts:
    if opt == "-b":
        print('Computing bandwidth.')
        compute_bandwidth=True
    elif opt == "-g":
        grid_file=arg
    elif opt == "-p":
        plots=True
    elif opt == "-s":
        sample=True
        sample_file=arg
        print('Computing samples and storing in file',sample_file,'.')
    elif opt == "-h":
        print('-b: compute bandwidth')
        print('-g <file>: specify tensor_grid filename as <file>')
        print('-h: show help')
        print('-p: make plots')
        print('-s <file>: compute samples and store in <file>')
    else:
        print('Incorrect option.')
        quit()

if plots==False:
        
    f=open("ligo_data","r" )
    lines=f.readlines()

    mass1=[]
    mass2=[]
    lambda1=[]
    lambda2=[]

    ###### making an array from given data ########
    for x in lines:
        x=x.split(' ')
        mass1.append(x[2])
        mass2.append(x[3])    
        lambda1.append(x[4])
        lambda2.append(x[5])    
    f.close()

    # Remove header row
    mass1.pop(0)
    mass2.pop(0)
    lambda1.pop(0)
    lambda2.pop(0)

    # Create numpy float arrays
    m1=np.asarray(mass1).astype(np.float)
    m2=np.asarray(mass2).astype(np.float)
    lam1=np.asarray(lambda1).astype(np.float)
    lam2=np.asarray(lambda2).astype(np.float)

    # Compute M_chirp, q, and lambda_tilde
    M_chirp=[]
    q=[]
    lambda_tilde=[]

    for i in range(0, len(m1)):
        M_chirp.append(((m1[i]*m2[i])**(3.0/5.0))/(m1[i]+m2[i])**(1.0/5.0))
        q.append(m2[i]/m1[i])
        eta=(m1[i]*m2[i])/((m1[i]+m2[i])*(m1[i]+m2[i]))
        Lt=(8.0/13.0)*((1+7*eta-31*eta*eta)*(lam1[i]+lam2[i])+
                       (1-4*eta)**(0.5)*(1+9*eta-11*eta*eta)*(lam1[i]-lam2[i]))
        lambda_tilde.append(Lt)

    # Compute mean and standard deviation

    M_chirp_mean=np.mean(M_chirp)
    M_chirp_std=np.std(M_chirp)

    q_mean=np.mean(q)
    q_std=np.std(q)

    lambda_tilde_mean=np.mean(lambda_tilde)
    lambda_tilde_std=np.std(lambda_tilde)

    # Standardize data

    M_chirp_transf=[]
    q_transf=[]
    lambda_tilde_transf=[]

    for i in range(0, len(M_chirp)):
        M_chirp_transf.append((M_chirp[i]-M_chirp_mean)/M_chirp_std)
        q_transf.append((q[i]-q_mean)/q_std)
        lambda_tilde_transf.append((lambda_tilde[i]-lambda_tilde_mean)/
                                   lambda_tilde_std)

    # Stack original and transformed data
    stack=np.vstack([M_chirp_transf,q_transf,lambda_tilde_transf])
    orig_data_stack=np.vstack([M_chirp,q,lambda_tilde])

    # Perform KDE, optimizing bandwidth if necessary
    if compute_bandwidth:
        grid=GridSearchCV(KernelDensity(),
                            {'bandwidth': np.linspace(0.1, 0.3, 20)},
                            verbose=2)
        grid.fit(stack.T)
        print('Best score:',grid.best_score_)
        print('Best bandwidth:',grid.best_estimator_.bandwidth)
        bw=grid.best_estimator_.bandwidth
        kde=grid.best_estimator_
    else:
        bw=0.2158
        kde=KernelDensity(bandwidth=bw).fit(stack.T)
    
    # pdf for samples in data stack 
    '''
    pdf=np.exp(kde.score_samples(stack.T))
    
    dist=np.concatenate((orig_data_stack.T, np.array([pdf]).T), axis=1)
    
    np.savetxt('orig_dist.txt', dist, fmt='%.6f', delimiter='\t', 
        header="chirp_mass\tmass_ratio\tlambda_tildea\tprob", comments='')
    '''
    
    # If requested, sample the distribution
    
    if sample:
        
        sample=kde.sample(n_samples=64000)
    
        # Invert the samples
        sample_M_chirp=[]
        sample_q=[]
        sample_lambda_tilde=[]
    
        for i in range(0, len(sample)):
            sample_M_chirp.append(sample[i][0]*M_chirp_std+M_chirp_mean)
            sample_q.append(sample[i][1]*q_std+q_mean)
            sample_lambda_tilde.append(sample[i][2]*lambda_tilde_std+
                                       lambda_tilde_mean)
    
        sample_pdf=np.exp(kde.score_samples(sample))
        sample_stack=np.vstack([sample_M_chirp,sample_q,sample_lambda_tilde])
        sample_dist=np.concatenate((sample_stack.T,np.array([sample_pdf]).T),
                                   axis=1)
        np.savetxt(sample_file,sample_dist,fmt='%.6f',delimiter=' ', 
                   header="chirp_mass mass_ratio lambda_tilde prob",
                   comments='')
    
    # Create the 3 grids in the original variables
    mint=min(M_chirp)
    maxt=max(M_chirp)
    low=mint-(maxt-mint)/10.0
    high=maxt+(maxt-mint)/10.0
    grid_mchirp=np.linspace(low,high,40)

    mint=min(q)
    maxt=max(q)
    low=mint-(maxt-mint)/10.0
    high=1.0
    grid_q=np.linspace(low,high,40)

    mint=min(lambda_tilde)
    maxt=max(lambda_tilde)
    low=0
    high=maxt+(maxt-mint)/10.0
    grid_lambda_tilde=np.linspace(low,high,40)

    # The np.vstack() function requires both to have same dimension.
    # Here, "grid_cube" array is initialized with zeroes with the
    # appropriate dimensions. This row is removed below.
    grid_cube=np.array([0,0,0])

    for x in grid_mchirp:
        for y in grid_q:
            for z in grid_lambda_tilde:
                temp=[x,y,z]
                grid_cube=np.vstack([grid_cube,[temp]])

    # Remove initialization array
    grid_cube=np.delete(grid_cube,0,0)

    grid_stack=np.vstack([grid_mchirp,grid_q,grid_lambda_tilde])

    # Now create the three grids in the transformed variables

    grid_trans_m=[]
    grid_trans_q=[]
    grid_trans_l=[]

    for i in range(0,len(grid_cube)):
        grid_trans_m.append((grid_cube[i][0]-M_chirp_mean)/M_chirp_std)
        grid_trans_q.append((grid_cube[i][1]-q_mean)/q_std)
        grid_trans_l.append((grid_cube[i][2]-lambda_tilde_mean)/
                            lambda_tilde_std)

    grid_sample_stack=np.vstack([grid_trans_m,grid_trans_q,grid_trans_l])

    # Compute the probabilities
    grid_log_pdf=kde.score_samples(grid_sample_stack.T)
    grid_pdf=np.exp(kde.score_samples(grid_sample_stack.T))

    # Write the probabilities to a file

    hf=h5py.File(grid_file,'w')

    # Create two tensor_grid objects, one for the probability and
    # one for the log of the probability
    
    for k in range(0,2):

        if k==0:
            group=hf.create_group('kde_prob')
        else:
            group=hf.create_group('kde_log_prob')
    
        data_dset=group.create_dataset('data',(64000,),dtype='d')
        if k==0:
            for i in range(0,64000):
                data_dset[i]=grid_pdf[i]
        else:
            for i in range(0,64000):
                data_dset[i]=grid_log_pdf[i]
        
        grid_dset=group.create_dataset('grid',(120,),dtype='d')
        for i in range(0,40):
            grid_dset[i]=grid_mchirp[i]
        for i in range(0,40):
            grid_dset[i+40]=grid_q[i]
        for i in range(0,40):
            grid_dset[i+80]=grid_lambda_tilde[i]
        grid_set_dset=group.create_dataset('grid_set',(1,),dtype='i')
        grid_set_dset[0]=1
    
        # h5py creates strings which are padded by NULL characters instead of
        # being terminated by them. In order to mimic a NULL terminated
        # string, we make the string one character larger
        o2scl_type_dset=group.create_dataset("o2scl_type", (1,), dtype="S12")
        o2scl_type_dset[0]=np.string_('tensor_grid')
    
        rank_dset=group.create_dataset('rank',(1,),dtype='i')
        rank_dset[0]=3
        size_dset=group.create_dataset('size',(3,),dtype='i')
        size_dset[0]=40
        size_dset[1]=40
        size_dset[2]=40

else:

    os.system('o2graph -read '+grid_file+
              ' kde_prob -rearrange \"sum(0),index(2),index(1)\" '+
              '-to-table3d 0 1 prob -den-plot prob '+
              '-xtitle \"$ \\tilde{\Lambda} $\" -ytitle \" $ q $ \" '+
              '-save lt_q.pdf')

