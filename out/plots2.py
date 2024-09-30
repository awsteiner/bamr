import numpy as np
import h5py
import matplotlib.pyplot as plt
import seaborn as sns
plt.rc('font', size=12.0)

nsrc = ['6304', '6397', 'M13', 'M28', 'M30', 'wCen', 'X7', '1810b', '1724b', '1702', '0030', '0740']
pstr = ['J0453p','J0453c','J1906p','J1906c','B1534p',  'B1534c','B1913p','B1913c',
      'B2127p','B2127c','J0737A','J0737B','J1756p',  'J1756c','J1807p','J1807c',
      'J1518p','J1518c','J1811p','J1811c','J1829p',  'J1829c','J2045', 'J2053',
      'J1713', 'B1855', 'J0751', 'J1141', 'J1738',   'J1614', 'J0348', 'J2222', 'J2234',
      'J1949', 'J1012', 'J0437', 'J1909', 'J1802',   'J1911', 'J2043', 'J0337', 'J1946',
      'J1918', 'J1600', 'J0621', 'B2303', 'J0024',   'J0514', 'B1516', 'J1748I','J1748J',
      'B1802', 'B1911', 'J0740', 'CygX2', 'XTEJ2123','4U1822','HerX1', '2S0921']

m_nsp, r_nsp, m_src, r_src = [], [], [], []
for k, star in enumerate(nsrc):
    m_src.append('Mns_'+nsrc[k])
for k, star in enumerate(nsrc):
    r_src.append('Rns_'+nsrc[k])
for k, star in enumerate(pstr):
    m_nsp.append('M_'+pstr[k])

fig, ax = plt.subplots(2, 2)
ax=ax.flatten()

model = ['ml', 'mp', 'nl', 'np']
for k in range(4):
    mchain = h5py.File(model[k]+'_all','r')['markov_chain_0']
    mmax = np.array(mchain['data']['M_max'])
    nrows = len(mmax)
    r = np.zeros((100, nrows))
    for i in range(100):
        r[i] = np.array(mchain['data']['R_'+str(i)])
    r=r.T
    m = np.arange(0.2, 0.2+((3-0.2)/99)*r.shape[1], ((3-0.2)/99))
    M_gw, R_gw = np.zeros((4, nrows)), np.zeros((4, nrows))
    M_gw[0] = np.array(mchain['data/m1_gw17'])
    M_gw[1] = np.array(mchain['data/m2_gw17'])
    M_gw[2] = np.array(mchain['data/m1_gw19'])
    M_gw[3] = np.array(mchain['data/m2_gw19'])
    M_nsp, R_nsp = np.zeros((len(m_nsp), nrows)), np.zeros((len(m_nsp), nrows))
    M_src, R_src = np.zeros((len(m_src), nrows)), np.zeros((len(m_src), nrows))

    for j in range(4):
        for i in range(nrows):
            R_gw[j][i] = np.interp(M_gw[j][i], m, r[i])
    for j in range(len(m_nsp)):
        M_nsp[j] = np.array(mchain['data/'+m_nsp[j]])
        for i in range(nrows):
            R_nsp[j][i] = np.interp(M_nsp[j][i], m, r[i])
    for j in range(len(m_src)):
        M_src[j] = np.array(mchain['data/'+m_src[j]])
        R_src[j] = np.array(mchain['data/'+r_src[j]])

    M = np.concatenate((M_gw, M_nsp, M_src), axis=0).flatten()
    R = np.concatenate((R_gw, R_nsp, R_src), axis=0).flatten()

    h = sns.kdeplot(x=R, y=M, ax=ax[k], fill=True, cmap='viridis')
    cb = fig.colorbar(h.collections[0], ax=ax[k])
    cb.locator = plt.MaxNLocator(nbins=5)
    cb.update_ticks()
    if k==0:
        ax[k].text(0.05, 0.95, 'ML', transform=ax[k].transAxes, \
                   verticalalignment='top')
    elif k==1:
        ax[k].text(0.05, 0.95, 'MP', transform=ax[k].transAxes, \
                   verticalalignment='top')
    elif k==2:
        ax[k].text(0.05, 0.95, 'NL', transform=ax[k].transAxes, \
                   verticalalignment='top')
    elif k==3:
        ax[k].text(0.05, 0.95, 'NP', transform=ax[k].transAxes, \
                   verticalalignment='top')
    ax[k].set_xlim(11, 13)
    ax[k].set_ylim(min(M), 2.15)
    ax[k].minorticks_on()
    ax[k].tick_params('both', length=10, width=1, which='major')
    ax[k].tick_params('both', length=5,  width=1, which='minor')
    ax[k].xaxis.set_major_locator(plt.MaxNLocator(nbins=3))
    ax[k].yaxis.set_major_locator(plt.MaxNLocator(nbins=3))

fig.text(0.5, 0.04, r'$R$ [km]', ha='center', va='center')
fig.text(0.04, 0.54, r'$M$ [$M_\odot$]', ha='center', va='center', \
         rotation='vertical')
fig.tight_layout()
plt.show()