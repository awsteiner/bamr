#! /usr/bin/bash

read -p 'Enter filename: [or path/to/file] > ' file
is_aff_inv=$(o2graph -read $file aff_inv -o)
echo -e "aff inv? $is_aff_inv"
read -p 'Plot by group: (1) EOS (2) LIGO (3) Population (4) P-E/M-R > ' plotgroup
if [ $plotgroup -eq 3 ]
then 
    read -p 'Population parameter type: (1) Distribution (2) True Mass > ' partype
    if [ $partype -eq 2 ]
    then 
        read -p 'Population type: (1) DNS (2) NSWD (3) HMXB (4) LMXB > ' poptype
    fi
fi
if [ $plotgroup -eq 4 ]
then 
    echo -e 'Making EOS and M-R plots...'
    o2graph -xtitle "Energy Density [fm$ ^{-4}$]" -ytitle "Pressure [fm$ ^{-4}$]" -plotv "grid:0.3,10.0,(10.0-0.3)/99" "hdf5:$file:markov_chain_0:*:P_*" marker=. -show &
    o2graph -xtitle "Radius [km]" -ytitle "M [$ M_\odot$]" -plotv "hdf5:$file:markov_chain_0:*:R_*" "grid:0.2,3.0,(3.0-0.2)/99" marker=. -show &
    echo 'Done.'
    exit
fi 
read -p 'Plot type: (1) Trace (2) Histogram (3) Likelihood > ' plottype 
if [ $plottype -eq 2 ]
then 
    read -p 'Enter bin size > ' binsz
fi

if [ $plotgroup -eq 1 ]
then 
    if [ $plottype -eq 1 ]
    then 
        echo -e 'Making trace plots of EOS parameters...'
        o2graph -read $file markov_chain_0 -ytitle a -xtitle steps -plot1 a -show &
        o2graph -read $file markov_chain_0 -ytitle alpha -xtitle steps -plot1 alpha -show &
        o2graph -read $file markov_chain_0 -ytitle param_S -xtitle steps -plot1 param_S -show &
        o2graph -read $file markov_chain_0 -ytitle param_L -xtitle steps -plot1 param_L -show &
        o2graph -read $file markov_chain_0 -ytitle exp1 -xtitle steps -plot1 exp1 -show &
        o2graph -read $file markov_chain_0 -ytitle trans1 -xtitle steps -plot1 trans1 -show &
        o2graph -read $file markov_chain_0 -ytitle exp2 -xtitle steps -plot1 exp2 -show &
        o2graph -read $file markov_chain_0 -ytitle trans2 -xtitle steps -plot1 trans2 -show &
        o2graph -read $file markov_chain_0 -ytitle exp3 -xtitle steps -plot1 exp3 -show &
    elif [ $plottype -eq 2 ]
    then 
        echo -e 'Making histogram plots of EOS parameters...'
        o2graph -read $file markov_chain_0 -xtitle a -ytitle frequency -to-hist a $binsz -plot -show &
        o2graph -read $file markov_chain_0 -xtitle alpha -ytitle frequency -to-hist alpha $binsz -plot -show &
        o2graph -read $file markov_chain_0 -xtitle param_S -ytitle frequency -to-hist param_S $binsz -plot -show &
        o2graph -read $file markov_chain_0 -xtitle param_L -ytitle frequency -to-hist param_L $binsz -plot -show &
        o2graph -read $file markov_chain_0 -xtitle exp1 -ytitle frequency -to-hist exp1 $binsz -plot -show &
        o2graph -read $file markov_chain_0 -xtitle trans1 -ytitle frequency -to-hist trans1 $binsz -plot -show &
        o2graph -read $file markov_chain_0 -xtitle exp2 -ytitle frequency -to-hist exp2 $binsz -plot -show &
        o2graph -read $file markov_chain_0 -xtitle trans2 -ytitle frequency -to-hist trans2 $binsz -plot -show &
        o2graph -read $file markov_chain_0 -xtitle exp3 -ytitle frequency -to-hist exp3 $binsz -plot -show &
    elif [ $plottype -eq 3 ]
    then 
        echo -e 'Making likelihood plots of EOS parameters...'
        o2graph -read $file markov_chain_0 -function 1 n -xtitle a -ytitle log_wgt -scatter a log_wgt n log_wgt -show &
        o2graph -read $file markov_chain_0 -function 1 n -xtitle alpha -ytitle log_wgt -scatter alpha log_wgt n log_wgt -show &
        o2graph -read $file markov_chain_0 -function 1 n -xtitle param_S -ytitle log_wgt -scatter param_S log_wgt n log_wgt -show &
        o2graph -read $file markov_chain_0 -function 1 n -xtitle param_L -ytitle log_wgt -scatter param_L log_wgt n log_wgt -show &
        o2graph -read $file markov_chain_0 -function 1 n -xtitle exp1 -ytitle log_wgt -scatter exp1 log_wgt n log_wgt -show &
        o2graph -read $file markov_chain_0 -function 1 n -xtitle trans1 -ytitle log_wgt -scatter trans1 log_wgt n log_wgt -show &
        o2graph -read $file markov_chain_0 -function 1 n -xtitle exp2 -ytitle log_wgt -scatter exp2 log_wgt n log_wgt -show &
        o2graph -read $file markov_chain_0 -function 1 n -xtitle trans2 -ytitle log_wgt -scatter trans2 log_wgt n log_wgt -show &
        o2graph -read $file markov_chain_0 -function 1 n -xtitle exp3 -ytitle log_wgt -scatter exp3 log_wgt n log_wgt -show &
    else 
        echo -e 'Invalid choice (exit).'
        exit
    fi 
    echo 'Done.'


elif [ $plotgroup -eq 2 ] 
then 
    if [ $plottype -eq 1 ]
    then 
        echo -e 'Making trace plots of LIGO parameters...'
        o2graph -read $file markov_chain_0 -ytitle M_chirp_det -xtitle steps -plot1 M_chirp_det -show &
        o2graph -read $file markov_chain_0 -ytitle q -xtitle steps -plot1 q -show &
        o2graph -read $file markov_chain_0 -ytitle z_cdf -xtitle steps -plot1 z_cdf -show &
    elif [ $plottype -eq 2 ]
    then 
        echo -e 'Making histogram plots of LIGO parameters...'
        o2graph -read $file markov_chain_0 -xtitle M_chirp_det -ytitle frequency -to-hist M_chirp_det $binsz -plot -show &
        o2graph -read $file markov_chain_0 -xtitle q -ytitle frequency -to-hist q $binsz -plot -show &
        o2graph -read $file markov_chain_0 -xtitle z_cdf -ytitle frequency -to-hist z_cdf $binsz -plot -show &
    elif [ $plottype -eq 3 ]
    then 
        echo -e 'Making likelihood plots of LIGO parameters...'
        o2graph -read $file markov_chain_0 -function 1 n -xtitle M_chirp_det -ytitle log_wgt -scatter M_chirp_det log_wgt n log_wgt -show &
        o2graph -read $file markov_chain_0 -function 1 n -xtitle q -ytitle log_wgt -scatter q log_wgt n log_wgt -show &
        o2graph -read $file markov_chain_0 -function 1 n -xtitle z_cdf -ytitle log_wgt -scatter z_cdf log_wgt n log_wgt -show &
    else 
        echo -e 'Invalid choice (exit).'
        exit
    fi 
    echo 'Done.'


elif [ $plotgroup -eq 3 ] && [ $partype -eq 1 ]
then 
    if [ $plottype -eq 1 ]
    then 
        echo -e 'Making trace plots of distribution parameters...'
        o2graph -read $file markov_chain_0 -ytitle mean_NS -xtitle steps -plot1 mean_NS -show &
        o2graph -read $file markov_chain_0 -ytitle log10_width_NS -xtitle steps -plot1 log10_width_NS -show &
        o2graph -read $file markov_chain_0 -ytitle skewness_NS -xtitle steps -plot1 skewness_NS -show &
        o2graph -read $file markov_chain_0 -ytitle mean_WD -xtitle steps -plot1 mean_WD -show &
        o2graph -read $file markov_chain_0 -ytitle log10_width_WD -xtitle steps -plot1 log10_width_WD -show &
        o2graph -read $file markov_chain_0 -ytitle skewness_WD -xtitle steps -plot1 skewness_WD -show &
        o2graph -read $file markov_chain_0 -ytitle mean_HMS -xtitle steps -plot1 mean_HMS -show &
        o2graph -read $file markov_chain_0 -ytitle log10_width_HMS -xtitle steps -plot1 log10_width_HMS -show &
        o2graph -read $file markov_chain_0 -ytitle skewness_HMS -xtitle steps -plot1 skewness_HMS -show &
        o2graph -read $file markov_chain_0 -ytitle mean_LMS -xtitle steps -plot1 mean_LMS -show &
        o2graph -read $file markov_chain_0 -ytitle log10_width_LMS -xtitle steps -plot1 log10_width_LMS -show &
        o2graph -read $file markov_chain_0 -ytitle skewness_LMS -xtitle steps -plot1 skewness_LMS -show &
    elif [ $plottype -eq 2 ]
    then 
        echo -e 'Making histogram plots of distribution parameters...'
        o2graph -read $file markov_chain_0 -xtitle mean_NS -ytitle frequency -to-hist mean_NS $binsz -plot -show &
        o2graph -read $file markov_chain_0 -xtitle log10_width_NS -ytitle frequency -to-hist log10_width_NS $binsz -plot -show &
        o2graph -read $file markov_chain_0 -xtitle skewness_NS -ytitle frequency -to-hist skewness_NS $binsz -plot -show &
        o2graph -read $file markov_chain_0 -xtitle mean_WD -ytitle frequency -to-hist mean_WD $binsz -plot -show &
        o2graph -read $file markov_chain_0 -xtitle log10_width_WD -ytitle frequency -to-hist log10_width_WD $binsz -plot -show &
        o2graph -read $file markov_chain_0 -xtitle skewness_WD -ytitle frequency -to-hist skewness_WD $binsz -plot -show &
        o2graph -read $file markov_chain_0 -xtitle mean_HMS -ytitle frequency -to-hist mean_HMS $binsz -plot -show &
        o2graph -read $file markov_chain_0 -xtitle log10_width_HMS -ytitle frequency -to-hist log10_width_HMS $binsz -plot -show &
        o2graph -read $file markov_chain_0 -xtitle skewness_HMS -ytitle frequency -to-hist skewness_HMS $binsz -plot -show &
        o2graph -read $file markov_chain_0 -xtitle mean_LMS -ytitle frequency -to-hist mean_LMS $binsz -plot -show &
        o2graph -read $file markov_chain_0 -xtitle log10_width_LMS -ytitle frequency -to-hist log10_width_LMS $binsz -plot -show &
        o2graph -read $file markov_chain_0 -xtitle skewness_LMS -ytitle frequency -to-hist skewness_LMS $binsz -plot -show &
    elif [ $plottype -eq 3 ]
    then 
        echo -e 'Making likelihood plots of distribution parameters...'
        o2graph -read $file markov_chain_0 -function 1 n -xtitle mean_NS -ytitle log_wgt -scatter mean_NS log_wgt n log_wgt -show &
        o2graph -read $file markov_chain_0 -function 1 n -xtitle log10_width_NS -ytitle log_wgt -scatter log10_width_NS log_wgt n log_wgt -show &
        o2graph -read $file markov_chain_0 -function 1 n -xtitle skewness_NS -ytitle log_wgt -scatter skewness_NS log_wgt n log_wgt -show &
        o2graph -read $file markov_chain_0 -function 1 n -xtitle mean_WD -ytitle log_wgt -scatter mean_WD log_wgt n log_wgt -show &
        o2graph -read $file markov_chain_0 -function 1 n -xtitle log10_width_WD -ytitle log_wgt -scatter log10_width_WD log_wgt n log_wgt -show &
        o2graph -read $file markov_chain_0 -function 1 n -xtitle skewness_WD -ytitle log_wgt -scatter skewness_WD log_wgt n log_wgt -show &
        o2graph -read $file markov_chain_0 -function 1 n -xtitle mean_HMS -ytitle log_wgt -scatter mean_HMS log_wgt n log_wgt -show &
        o2graph -read $file markov_chain_0 -function 1 n -xtitle log10_width_HMS -ytitle log_wgt -scatter log10_width_HMS log_wgt n log_wgt -show &
        o2graph -read $file markov_chain_0 -function 1 n -xtitle skewness_HMS -ytitle log_wgt -scatter skewness_HMS log_wgt n log_wgt -show &
        o2graph -read $file markov_chain_0 -function 1 n -xtitle mean_LMS -ytitle log_wgt -scatter mean_LMS log_wgt n log_wgt -show &
        o2graph -read $file markov_chain_0 -function 1 n -xtitle log10_width_LMS -ytitle log_wgt -scatter log10_width_LMS log_wgt n log_wgt -show &
        o2graph -read $file markov_chain_0 -function 1 n -xtitle skewness_LMS -ytitle log_wgt -scatter skewness_LMS log_wgt n log_wgt -show &
    else 
        echo -e 'Invalid choice (exit).'
        exit
    fi 
    echo 'Done.'


elif [ $plotgroup -eq 3 ] && [ $partype -eq 2 ]
then 
    if [ $plottype -eq 1 ]
    then 
        if [ $poptype -eq 1 ]
        then 
            echo -e 'Making trace plots of DNS true-mass parameters...'
            o2graph -read $file markov_chain_0 -ytitle M_J0453p -xtitle steps -plot1 M_J0453p -show &
            o2graph -read $file markov_chain_0 -ytitle M_J0453c -xtitle steps -plot1 M_J0453c -show &
            o2graph -read $file markov_chain_0 -ytitle M_J1906p -xtitle steps -plot1 M_J1906p -show &
            o2graph -read $file markov_chain_0 -ytitle M_J1906c -xtitle steps -plot1 M_J1906c -show &
            o2graph -read $file markov_chain_0 -ytitle M_B1534p -xtitle steps -plot1 M_B1534p -show &
            o2graph -read $file markov_chain_0 -ytitle M_B1534c -xtitle steps -plot1 M_B1534c -show &
            o2graph -read $file markov_chain_0 -ytitle M_B1913p -xtitle steps -plot1 M_B1913p -show &
            o2graph -read $file markov_chain_0 -ytitle M_B1913c -xtitle steps -plot1 M_B1913c -show &
            o2graph -read $file markov_chain_0 -ytitle M_B2127p -xtitle steps -plot1 M_B2127p -show &
            o2graph -read $file markov_chain_0 -ytitle M_B2127c -xtitle steps -plot1 M_B2127c -show &
            o2graph -read $file markov_chain_0 -ytitle M_J0737A -xtitle steps -plot1 M_J0737A -show &
            o2graph -read $file markov_chain_0 -ytitle M_J0737B -xtitle steps -plot1 M_J0737B -show &
            o2graph -read $file markov_chain_0 -ytitle M_J1756p -xtitle steps -plot1 M_J1756p -show &
            o2graph -read $file markov_chain_0 -ytitle M_J1756c -xtitle steps -plot1 M_J1756c -show &
            o2graph -read $file markov_chain_0 -ytitle M_J1807p -xtitle steps -plot1 M_J1807p -show &
            o2graph -read $file markov_chain_0 -ytitle M_J1807c -xtitle steps -plot1 M_J1807c -show &
            o2graph -read $file markov_chain_0 -ytitle M_J1518p -xtitle steps -plot1 M_J1518p -show &
            o2graph -read $file markov_chain_0 -ytitle M_J1518c -xtitle steps -plot1 M_J1518c -show &
            o2graph -read $file markov_chain_0 -ytitle M_J1811p -xtitle steps -plot1 M_J1811p -show &
            o2graph -read $file markov_chain_0 -ytitle M_J1811c -xtitle steps -plot1 M_J1811c -show &
            o2graph -read $file markov_chain_0 -ytitle M_J1829p -xtitle steps -plot1 M_J1829p -show &
            o2graph -read $file markov_chain_0 -ytitle M_J1829c -xtitle steps -plot1 M_J1829c -show &
        elif [ $poptype -eq 2 ]
        then 
            echo -e 'Making trace plots of NSWD true-mass parameters...'
            o2graph -read $file markov_chain_0 -ytitle M_J2045 -xtitle steps -plot1 M_J2045 -show &
            o2graph -read $file markov_chain_0 -ytitle M_J2053 -xtitle steps -plot1 M_J2053 -show &
            o2graph -read $file markov_chain_0 -ytitle M_J1713 -xtitle steps -plot1 M_J1713 -show &
            o2graph -read $file markov_chain_0 -ytitle M_B1855 -xtitle steps -plot1 M_B1855 -show &
            o2graph -read $file markov_chain_0 -ytitle M_J0751 -xtitle steps -plot1 M_J0751 -show &
            o2graph -read $file markov_chain_0 -ytitle M_J1141 -xtitle steps -plot1 M_J1141 -show &
            o2graph -read $file markov_chain_0 -ytitle M_J1738 -xtitle steps -plot1 M_J1738 -show &
            o2graph -read $file markov_chain_0 -ytitle M_J1614 -xtitle steps -plot1 M_J1614 -show &
            o2graph -read $file markov_chain_0 -ytitle M_J0348 -xtitle steps -plot1 M_J0348 -show &
            o2graph -read $file markov_chain_0 -ytitle M_J2222 -xtitle steps -plot1 M_J2222 -show &
            o2graph -read $file markov_chain_0 -ytitle M_J2234 -xtitle steps -plot1 M_J2234 -show &
            o2graph -read $file markov_chain_0 -ytitle M_J1949 -xtitle steps -plot1 M_J1949 -show &
            o2graph -read $file markov_chain_0 -ytitle M_J1012 -xtitle steps -plot1 M_J1012 -show &
            o2graph -read $file markov_chain_0 -ytitle M_J0437 -xtitle steps -plot1 M_J0437 -show &
            o2graph -read $file markov_chain_0 -ytitle M_J1909 -xtitle steps -plot1 M_J1909 -show &
            o2graph -read $file markov_chain_0 -ytitle M_J1802 -xtitle steps -plot1 M_J1802 -show &
            o2graph -read $file markov_chain_0 -ytitle M_J1911 -xtitle steps -plot1 M_J1911 -show &
            o2graph -read $file markov_chain_0 -ytitle M_J2043 -xtitle steps -plot1 M_J2043 -show &
            o2graph -read $file markov_chain_0 -ytitle M_J0337 -xtitle steps -plot1 M_J0337 -show &
            o2graph -read $file markov_chain_0 -ytitle M_J1946 -xtitle steps -plot1 M_J1946 -show &
            o2graph -read $file markov_chain_0 -ytitle M_J1918 -xtitle steps -plot1 M_J1918 -show &
            o2graph -read $file markov_chain_0 -ytitle M_J1600 -xtitle steps -plot1 M_J1600 -show &
            o2graph -read $file markov_chain_0 -ytitle M_J0621 -xtitle steps -plot1 M_J0621 -show &
            o2graph -read $file markov_chain_0 -ytitle M_B2303 -xtitle steps -plot1 M_B2303 -show &
            o2graph -read $file markov_chain_0 -ytitle M_J0024 -xtitle steps -plot1 M_J0024 -show &
            o2graph -read $file markov_chain_0 -ytitle M_J0514 -xtitle steps -plot1 M_J0514 -show &
            o2graph -read $file markov_chain_0 -ytitle M_B1516 -xtitle steps -plot1 M_B1516 -show &
            o2graph -read $file markov_chain_0 -ytitle M_J1748I -xtitle steps -plot1 M_J1748I -show &
            o2graph -read $file markov_chain_0 -ytitle M_J1748J -xtitle steps -plot1 M_J1748J -show &
            o2graph -read $file markov_chain_0 -ytitle M_B1802 -xtitle steps -plot1 M_B1802 -show &
            o2graph -read $file markov_chain_0 -ytitle M_B1911 -xtitle steps -plot1 M_B1911 -show &
            o2graph -read $file markov_chain_0 -ytitle M_J0740 -xtitle steps -plot1 M_J0740 -show &
        elif [ $poptype -eq 3 ]
        then 
            echo -e 'Making trace plots of HMXB true-mass parameters...'
            o2graph -read $file markov_chain_0 -ytitle M_4U1700 -xtitle steps -plot1 M_4U1700 -show &
            o2graph -read $file markov_chain_0 -ytitle M_SMCX1 -xtitle steps -plot1 M_SMCX1 -show &
            o2graph -read $file markov_chain_0 -ytitle M_CenX3 -xtitle steps -plot1 M_CenX3 -show &
            o2graph -read $file markov_chain_0 -ytitle M_OAO1657 -xtitle steps -plot1 M_OAO1657 -show &
            o2graph -read $file markov_chain_0 -ytitle M_J013236 -xtitle steps -plot1 M_J013236 -show &
            o2graph -read $file markov_chain_0 -ytitle M_VelaX1 -xtitle steps -plot1 M_VelaX1 -show &
            o2graph -read $file markov_chain_0 -ytitle M_4U1538 -xtitle steps -plot1 M_4U1538 -show &
            o2graph -read $file markov_chain_0 -ytitle M_LMCX4 -xtitle steps -plot1 M_LMCX4 -show &
            o2graph -read $file markov_chain_0 -ytitle M_EXO1722 -xtitle steps -plot1 M_EXO1722 -show &
            o2graph -read $file markov_chain_0 -ytitle M_SAXJ1802 -xtitle steps -plot1 M_SAXJ1802 -show &
            o2graph -read $file markov_chain_0 -ytitle M_XTEJ1855 -xtitle steps -plot1 M_XTEJ1855 -show &
        elif [ $poptype -eq 4 ]
        then 
            echo -e 'Making trace plots of LMXB true-mass parameters...'
            o2graph -read $file markov_chain_0 -ytitle M_CygX2 -xtitle steps -plot1 M_CygX2 -show &
            o2graph -read $file markov_chain_0 -ytitle M_XTEJ2123 -xtitle steps -plot1 M_XTEJ2123 -show &
            o2graph -read $file markov_chain_0 -ytitle M_4U1822 -xtitle steps -plot1 M_4U1822 -show &
            o2graph -read $file markov_chain_0 -ytitle M_HerX1 -xtitle steps -plot1 M_HerX1 -show &
            o2graph -read $file markov_chain_0 -ytitle M_2S0921 -xtitle steps -plot1 M_2S0921 -show &
        else 
            echo -e 'Invalid choice (exit).'
            exit
        fi 
        echo 'Done.'
    
    elif [ $plottype -eq 2 ]
    then 
        if [ $poptype -eq 1 ]
        then 
            echo -e 'Making histogram plots of DNS true-mass parameters...'
            o2graph -read $file markov_chain_0 -xtitle M_J0453p -ytitle frequency -to-hist M_J0453p $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_J0453c -ytitle frequency -to-hist M_J0453c $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_J1906p -ytitle frequency -to-hist M_J1906p $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_J1906c -ytitle frequency -to-hist M_J1906c $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_B1534p -ytitle frequency -to-hist M_B1534p $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_B1534c -ytitle frequency -to-hist M_B1534c $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_B1913p -ytitle frequency -to-hist M_B1913p $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_B1913c -ytitle frequency -to-hist M_B1913c $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_B2127p -ytitle frequency -to-hist M_B2127p $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_B2127c -ytitle frequency -to-hist M_B2127c $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_J0737A -ytitle frequency -to-hist M_J0737A $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_J0737B -ytitle frequency -to-hist M_J0737B $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_J1756p -ytitle frequency -to-hist M_J1756p $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_J1756c -ytitle frequency -to-hist M_J1756c $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_J1807p -ytitle frequency -to-hist M_J1807p $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_J1807c -ytitle frequency -to-hist M_J1807c $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_J1518p -ytitle frequency -to-hist M_J1518p $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_J1518c -ytitle frequency -to-hist M_J1518c $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_J1811p -ytitle frequency -to-hist M_J1811p $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_J1811c -ytitle frequency -to-hist M_J1811c $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_J1829p -ytitle frequency -to-hist M_J1829p $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_J1829c -ytitle frequency -to-hist M_J1829c $binsz -plot -show &
        elif [ $poptype -eq 2 ]
        then 
            echo -e 'Making histogram plots of NSWD true-mass parameters...'
            o2graph -read $file markov_chain_0 -xtitle M_J2045 -ytitle frequency -to-hist M_J2045 $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_J2053 -ytitle frequency -to-hist M_J2053 $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_J1713 -ytitle frequency -to-hist M_J1713 $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_B1855 -ytitle frequency -to-hist M_B1855 $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_J0751 -ytitle frequency -to-hist M_J0751 $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_J1141 -ytitle frequency -to-hist M_J1141 $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_J1738 -ytitle frequency -to-hist M_J1738 $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_J1614 -ytitle frequency -to-hist M_J1614 $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_J0348 -ytitle frequency -to-hist M_J0348 $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_J2222 -ytitle frequency -to-hist M_J2222 $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_J2234 -ytitle frequency -to-hist M_J2234 $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_J1949 -ytitle frequency -to-hist M_J1949 $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_J1012 -ytitle frequency -to-hist M_J1012 $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_J0437 -ytitle frequency -to-hist M_J0437 $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_J1909 -ytitle frequency -to-hist M_J1909 $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_J1802 -ytitle frequency -to-hist M_J1802 $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_J1911 -ytitle frequency -to-hist M_J1911 $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_J2043 -ytitle frequency -to-hist M_J2043 $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_J0337 -ytitle frequency -to-hist M_J0337 $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_J1946 -ytitle frequency -to-hist M_J1946 $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_J1918 -ytitle frequency -to-hist M_J1918 $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_J1600 -ytitle frequency -to-hist M_J1600 $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_J0621 -ytitle frequency -to-hist M_J0621 $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_B2303 -ytitle frequency -to-hist M_B2303 $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_J0024 -ytitle frequency -to-hist M_J0024 $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_J0514 -ytitle frequency -to-hist M_J0514 $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_B1516 -ytitle frequency -to-hist M_B1516 $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_J1748I -ytitle frequency -to-hist M_J1748I $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_J1748J -ytitle frequency -to-hist M_J1748J $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_B1802 -ytitle frequency -to-hist M_B1802 $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_B1911 -ytitle frequency -to-hist M_B1911 $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_J0740 -ytitle frequency -to-hist M_J0740 $binsz -plot -show &
        elif [ $poptype -eq 3 ]
        then 
            echo -e 'Making histogram plots of HMXB true-mass parameters...'
            o2graph -read $file markov_chain_0 -xtitle M_4U1700 -ytitle frequency -to-hist M_4U1700 $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_SMCX1 -ytitle frequency -to-hist M_SMCX1 $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_CenX3 -ytitle frequency -to-hist M_CenX3 $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_OAO1657 -ytitle frequency -to-hist M_OAO1657 $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_J013236 -ytitle frequency -to-hist M_J013236 $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_VelaX1 -ytitle frequency -to-hist M_VelaX1 $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_4U1538 -ytitle frequency -to-hist M_4U1538 $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_LMCX4 -ytitle frequency -to-hist M_LMCX4 $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_EXO1722 -ytitle frequency -to-hist M_EXO1722 $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_SAXJ1802 -ytitle frequency -to-hist M_SAXJ1802 $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_XTEJ1855 -ytitle frequency -to-hist M_XTEJ1855 $binsz -plot -show &
        elif [ $poptype -eq 4 ]
        then 
            echo -e 'Making histogram plots of LMXB true-mass parameters...'
            o2graph -read $file markov_chain_0 -xtitle M_CygX2 -ytitle frequency -to-hist M_CygX2 $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_XTEJ2123 -ytitle frequency -to-hist M_XTEJ2123 $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_4U1822 -ytitle frequency -to-hist M_4U1822 $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_HerX1 -ytitle frequency -to-hist M_HerX1 $binsz -plot -show &
            o2graph -read $file markov_chain_0 -xtitle M_2S0921 -ytitle frequency -to-hist M_2S0921 $binsz -plot -show &
        else 
            echo -e 'Invalid choice (exit).'
            exit
        fi 
        echo 'Done.'


    elif [ $plottype -eq 3 ]
    then 
        if [ $poptype -eq 1 ]
        then 
            echo -e 'Making likelihood plots of DNS true-mass parameters...'
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_J0453p -ytitle log_wgt -scatter M_J0453p log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_J0453c -ytitle log_wgt -scatter M_J0453c log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_J1906p -ytitle log_wgt -scatter M_J1906p log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_J1906c -ytitle log_wgt -scatter M_J1906c log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_B1534p -ytitle log_wgt -scatter M_B1534p log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_B1534c -ytitle log_wgt -scatter M_B1534c log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_B1913p -ytitle log_wgt -scatter M_B1913p log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_B1913c -ytitle log_wgt -scatter M_B1913c log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_B2127p -ytitle log_wgt -scatter M_B2127p log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_B2127c -ytitle log_wgt -scatter M_B2127c log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_J0737A -ytitle log_wgt -scatter M_J0737A log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_J0737B -ytitle log_wgt -scatter M_J0737B log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_J1756p -ytitle log_wgt -scatter M_J1756p log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_J1756c -ytitle log_wgt -scatter M_J1756c log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_J1807p -ytitle log_wgt -scatter M_J1807p log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_J1807c -ytitle log_wgt -scatter M_J1807c log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_J1518p -ytitle log_wgt -scatter M_J1518p log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_J1518c -ytitle log_wgt -scatter M_J1518c log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_J1811p -ytitle log_wgt -scatter M_J1811p log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_J1811c -ytitle log_wgt -scatter M_J1811c log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_J1829p -ytitle log_wgt -scatter M_J1829p log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_J1829c -ytitle log_wgt -scatter M_J1829c log_wgt n log_wgt -show &
        elif [ $poptype -eq 2 ]
        then 
            echo -e 'Making likelihood plots of NSWD true-mass parameters...'
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_J2045 -ytitle log_wgt -scatter M_J2045 log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_J2053 -ytitle log_wgt -scatter M_J2053 log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_J1713 -ytitle log_wgt -scatter M_J1713 log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_B1855 -ytitle log_wgt -scatter M_B1855 log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_J0751 -ytitle log_wgt -scatter M_J0751 log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_J1141 -ytitle log_wgt -scatter M_J1141 log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_J1738 -ytitle log_wgt -scatter M_J1738 log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_J1614 -ytitle log_wgt -scatter M_J1614 log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_J0348 -ytitle log_wgt -scatter M_J0348 log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_J2222 -ytitle log_wgt -scatter M_J2222 log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_J2234 -ytitle log_wgt -scatter M_J2234 log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_J1949 -ytitle log_wgt -scatter M_J1949 log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_J1012 -ytitle log_wgt -scatter M_J1012 log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_J0437 -ytitle log_wgt -scatter M_J0437 log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_J1909 -ytitle log_wgt -scatter M_J1909 log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_J1802 -ytitle log_wgt -scatter M_J1802 log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_J1911 -ytitle log_wgt -scatter M_J1911 log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_J2043 -ytitle log_wgt -scatter M_J2043 log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_J0337 -ytitle log_wgt -scatter M_J0337 log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_J1946 -ytitle log_wgt -scatter M_J1946 log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_J1918 -ytitle log_wgt -scatter M_J1918 log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_J1600 -ytitle log_wgt -scatter M_J1600 log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_J0621 -ytitle log_wgt -scatter M_J0621 log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_B2303 -ytitle log_wgt -scatter M_B2303 log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_J0024 -ytitle log_wgt -scatter M_J0024 log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_J0514 -ytitle log_wgt -scatter M_J0514 log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_B1516 -ytitle log_wgt -scatter M_B1516 log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_J1748I -ytitle log_wgt -scatter M_J1748I log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_J1748J -ytitle log_wgt -scatter M_J1748J log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_B1802 -ytitle log_wgt -scatter M_B1802 log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_B1911 -ytitle log_wgt -scatter M_B1911 log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_J0740 -ytitle log_wgt -scatter M_J0740 log_wgt n log_wgt -show &
        elif [ $poptype -eq 3 ]
        then 
            echo -e 'Making likelihood plots of HMXB true-mass parameters...'
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_4U1700 -ytitle log_wgt -scatter M_4U1700 log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_SMCX1 -ytitle log_wgt -scatter M_SMCX1 log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_CenX3 -ytitle log_wgt -scatter M_CenX3 log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_OAO1657 -ytitle log_wgt -scatter M_OAO1657 log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_J013236 -ytitle log_wgt -scatter M_J013236 log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_VelaX1 -ytitle log_wgt -scatter M_VelaX1 log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_4U1538 -ytitle log_wgt -scatter M_4U1538 log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_LMCX4 -ytitle log_wgt -scatter M_LMCX4 log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_EXO1722 -ytitle log_wgt -scatter M_EXO1722 log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_SAXJ1802 -ytitle log_wgt -scatter M_SAXJ1802 log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_XTEJ1855 -ytitle log_wgt -scatter M_XTEJ1855 log_wgt n log_wgt -show &
        elif [ $poptype -eq 4 ]
        then 
            echo -e 'Making likelihood plots of LMXB true-mass parameters...'
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_CygX2 -ytitle log_wgt -scatter M_CygX2 log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_XTEJ2123 -ytitle log_wgt -scatter M_XTEJ2123 log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_4U1822 -ytitle log_wgt -scatter M_4U1822 log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_HerX1 -ytitle log_wgt -scatter M_HerX1 log_wgt n log_wgt -show &
            o2graph -read $file markov_chain_0 -function 1 n -xtitle M_2S0921 -ytitle log_wgt -scatter M_2S0921 log_wgt n log_wgt -show &
        else 
            echo -e 'Invalid choice (exit).'
            exit
        fi 
        echo 'Done.'
    else 
        echo -e 'Invalid choice (exit).'
        exit
    fi 
else 
    echo -e 'Invalid choice (exit).'
    exit
fi 
