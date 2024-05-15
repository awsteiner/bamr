
model =
title =

INDEX = $(shell seq 0 99)
THREAD = $(shell seq 0 7)
WALKER = $(shell seq 0 63)
IWALK = $(shell seq 0 511)

new-file:
	acol -read mp2_31 markov_chain_0 $(foreach i, $(IWALK),-cat mp2_31) -internal mp3_31

# Base file - create either by concatenating, or copying, or processing raw files
file-main:
	acol -read out_$(model) markov_chain_0 -select-rows "mult>0.5" \
	-internal files/tmp_$(model).o2
	acol -h5-copy files/tmp_$(model).o2 files/main_$(model).o2
	rm -f files/tmp_$(model).o2

# Requires main_$(model).o2
file-esmr:
	acol -read files/main_$(model).o2 markov_chain_0 \
	$(foreach i,$(INDEX),-function P_$i*197.32698 Phc_$i) \
	-set obj_name ed_pr -select ":Phc_*" -internal files/esmr_$(model).o2 \
	-read files/esmr_$(model).o2 ed_pr -ser-hist-t3d \
	"grid:0.3*197.32698,10.0*197.32698,((10.0-0.3)*197.32698)/99.0" x ed \
	"grid:1.0,2.0e3,(2.0e3)^(1.0/99.0),log" auto pr \
	"Phc_*" prob -type -list \
	-set obj_name ed_pr_3d -internal files/esmr_$(model).o2 \
	-to-hist-2d prob \
	-contours frac 0.68 files/esmr_$(model).o2 c68_es \
	-contours frac 0.95 files/esmr_$(model).o2 c95_es \
	-read files/main_$(model).o2 markov_chain_0 -delete-col R_max \
	-set obj_name gm_r -select ":R_*" -internal files/esmr_$(model).o2 \
	-read files/esmr_$(model).o2 gm_r -ser-hist-t3d \
	"grid:0.2,3.0,(3.0-0.2)/99.0" y gm \
	"grid:9.0,15.0,(15.0-9.0)/99.0" auto r \
	"R_*" prob -type -list \
	-set obj_name gm_r_3d -internal files/esmr_$(model).o2 \
	-to-hist-2d prob \
	-contours frac 0.68 files/esmr_$(model).o2 c68_mr \
	-contours frac 0.95 files/esmr_$(model).o2 c95_mr

# Requires main_$(model).o2
file-main2:
	acol -read files/main_$(model).o2 markov_chain_0 \
	-function "10^log10_width_NS" width_NS \
	-function "10^log10_width_WD" width_WD \
	-function "10^log10_width_LMS" width_LM \
	-function "mean_LMS" mean_LM -delete-col mean_LMS \
	-function "skewness_LMS" skewness_LM \
	-delete-col skewness_LMS \
	-function "1/sqrt(2*3.1416)/width_NS" cf \
	-function "exp(-0.5*((1.4-mean_NS)/width_NS)^2)" pdf \
	-function "1+erf((1.4-mean_NS)*(skewness_NS)/width_NS/sqrt(2))" cdf \
	-function "cf*pdf*cdf" SN_NS \
	-function "1/sqrt(2*3.1416)/width_WD" cf \
	-function "exp(-0.5*((1.6-mean_WD)/width_WD)^2)" pdf \
	-function "1+erf((1.6-mean_WD)*(skewness_WD)/width_WD/sqrt(2))" cdf \
	-function "cf*pdf*cdf" SN_WD \
	-function "1/sqrt(2*3.1416)/width_LM" cf \
	-function "exp(-0.5*((1.7-mean_LM)/width_LM)^2)" pdf \
	-function "1+erf((1.7-mean_LM)*(skewness_LM)/width_LM/sqrt(2))" cdf \
	-function "cf*pdf*cdf" SN_LM \
	-function "log(dpdM)" log_dpdM \
	-function "P_23*197.32698" Phc_23 \
	-select ":mean_*" ":width_*" ":skewness_*" ":SN_*" \
	Phc_23 R_43 R_50 R_53 M_max log_dpdM log_wgt \
	-internal files/main2_$(model).o2

XLO = 0.0
XHI = 3.0
YLO = 0.0
YHI = 4.0
DX = 99.0

# Requires main2_$(model).o2
file-dist:
	acol -read files/main2_$(model).o2 markov_chain_0 \
	$(foreach i,$(INDEX),-function "1/sqrt(2*3.1416)/width_NS" cf \
	-function "exp(-0.5*(($i*0.03-mean_NS)/width_NS)^2)" pdf \
	-function "1+erf(($i*0.03-mean_NS)*(skewness_NS)/width_NS/sqrt(2))" cdf \
	-function "if($i*0.03<=M_max,cf*pdf*cdf,0.0)" SN_NS_$i) \
	$(foreach i,$(INDEX),-function "1/sqrt(2*3.1416)/width_WD" cf \
	-function "exp(-0.5*(($i*0.03-mean_WD)/width_WD)^2)" pdf \
	-function "1+erf(($i*0.03-mean_WD)*(skewness_WD)/width_WD/sqrt(2))" cdf \
	-function "if($i*0.03<=M_max,cf*pdf*cdf,0.0)" SN_WD_$i) \
	$(foreach i,$(INDEX),-function "1/sqrt(2*3.1416)/width_LM" cf \
	-function "exp(-0.5*(($i*0.03-mean_LM)/width_LM)^2)" pdf \
	-function "1+erf(($i*0.03-mean_LM)*(skewness_LM)/width_LM/sqrt(2))" cdf \
	-function "if($i*0.03<=M_max,cf*pdf*cdf,0.0)" SN_LM_$i) \
	-delete-col SN_NS -delete-col SN_WD -delete-col SN_LM \
	-select ":SN_*" -set obj_name gm_sn -internal files/dist_$(model).o2
	acol -read files/dist_$(model).o2 gm_sn -ser-hist-t3d \
	"grid:$(XLO),$(XHI),($(XHI)-$(XLO))/$(DX)" x gm \
	"grid:$(YLO),$(YHI),($(YHI)-$(YLO))/$(DX)" auto sn \
	"SN_NS_*" prob_ns -type -list \
	-set obj_name sn_ns -internal files/dist_$(model).o2 \
	-to-hist-2d prob_ns \
	-contours frac 0.68 files/dist_$(model).o2 c68_ns \
	-contours frac 0.95 files/dist_$(model).o2 c95_ns \
	-read files/dist_$(model).o2 gm_sn -ser-hist-t3d \
	"grid:$(XLO),$(XHI),($(XHI)-$(XLO))/$(DX)" x gm \
	"grid:$(YLO),$(YHI),($(YHI)-$(YLO))/$(DX)" auto sn \
	"SN_WD_*" prob_wd -type -list \
	-set obj_name sn_wd -internal files/dist_$(model).o2 \
	-to-hist-2d prob_wd \
	-contours frac 0.68 files/dist_$(model).o2 c68_wd \
	-contours frac 0.95 files/dist_$(model).o2 c95_wd \
	-read files/dist_$(model).o2 gm_sn -ser-hist-t3d \
	"grid:$(XLO),$(XHI),($(XHI)-$(XLO))/$(DX)" x gm \
	"grid:$(YLO),$(YHI),($(YHI)-$(YLO))/$(DX)" auto sn \
	"SN_LM_*" prob_lm -type -list \
	-set obj_name sn_lm -internal files/dist_$(model).o2 \
	-to-hist-2d prob_lm \
	-contours frac 0.68 files/dist_$(model).o2 c68_lm \
	-contours frac 0.95 files/dist_$(model).o2 c95_lm 

plot-dist:
	o2graph -read files/dist_$(model).o2 sn_ns -set colbar 1 \
	-den-plot prob_ns cmap=pink_r \
	-read files/dist_$(model).o2 c68_ns -plot lw=1,color=purple \
	-read files/dist_$(model).o2 c95_ns -plot lw=1,color=orange \
	-xtitle "Mass, $$ M $$ [$$ M_\odot $$]" \
	-ytitle "Mass Distribution, $$ \mathcal{P} $$" \
	-ttext 0.5 0.96 "$(model): DNS" \
	-save figs/sn_ns_$(model).png -show & 
	o2graph -read files/dist_$(model).o2 sn_wd -set colbar 1 \
	-den-plot prob_wd cmap=pink_r \
	-read files/dist_$(model).o2 c68_wd -plot lw=1,color=purple \
	-read files/dist_$(model).o2 c95_wd -plot lw=1,color=orange \
	-xtitle "Mass, $$ M $$ [$$ M_\odot $$]" \
	-ytitle "Mass Distribution, $$ \mathcal{P} $$" \
	-ttext 0.5 0.96 "$(model): NS-WD" \
	-save figs/sn_wd_$(model).png -show & 
	o2graph -read files/dist_$(model).o2 sn_lm -set colbar 1 \
	-den-plot prob_lm cmap=pink_r \
	-read files/dist_$(model).o2 c68_lm -plot lw=1,color=purple \
	-read files/dist_$(model).o2 c95_lm -plot lw=1,color=orange \
	-xtitle "Mass, $$ M $$ [$$ M_\odot $$]" \
	-ytitle "Mass Distribution, $$ \mathcal{P} $$" \
	-ttext 0.5 0.96 "$(model): LMXB" \
	-save figs/sn_lm_$(model).png -show & 

splot-eos1:
	o2graph -subplots 2 1 sharex=True \
	-subadj hspace=0.0,left=0.14,right=0.84,top=0.98,bottom=0.12 \
	-selax 0 -set logy 1 -read files/esmr_NL.o2 ed_pr_3d \
	-den-plot prob pcm=True,cmap=pink_r \
	-read files/esmr_NL.o2 c68_es -plot lw=1,color=purple \
	-read files/esmr_NL.o2 c95_es -plot lw=1,color=orange \
	-ttext 0.1 0.92 "NL" \
	-selax 1 -set logy 1 -read files/esmr_NP.o2 ed_pr_3d \
	-den-plot prob pcm=True,cmap=pink_r \
	-read files/esmr_NP.o2 c68_es -plot lw=1,color=purple \
	-read files/esmr_NP.o2 c95_es -plot lw=1,color=orange \
	-xtitle "Energy Density, $$ \epsilon~[\mathrm{MeV}/\mathrm{fm}^3] $$" \
	-ytitle ".$$ \qquad \qquad \qquad \qquad \qquad $$ Pressure, $$ P~[\mathrm{MeV}/\mathrm{fm}^3] $$" \
	-ttext 0.1 0.92 "NP" \
	-addcbar 0.86 0.12 0.03 0.86 -selax \
	-save figs/eos_s1.png -show &

splot-eos2:
	o2graph -subplots 2 1 sharex=True \
	-subadj hspace=0.0,left=0.14,right=0.84,top=0.98,bottom=0.12 \
	-selax 0 -set logy 1 -read files/esmr_ML.o2 ed_pr_3d \
	-den-plot prob pcm=True,cmap=pink_r \
	-read files/esmr_ML.o2 c68_es -plot lw=1,color=purple \
	-read files/esmr_ML.o2 c95_es -plot lw=1,color=orange \
	-ttext 0.1 0.92 "ML" \
	-selax 1 -set logy 1 -read files/esmr_MP.o2 ed_pr_3d \
	-den-plot prob pcm=True,cmap=pink_r \
	-read files/esmr_MP.o2 c68_es -plot lw=1,color=purple \
	-read files/esmr_MP.o2 c95_es -plot lw=1,color=orange \
	-xtitle "Energy Density, $$ \epsilon~[\mathrm{MeV}/\mathrm{fm}^3] $$" \
	-ytitle ".$$ \qquad \qquad \qquad \qquad \qquad $$ Pressure, $$ P~[\mathrm{MeV}/\mathrm{fm}^3] $$" \
	-ttext 0.1 0.92 "MP" \
	-addcbar 0.86 0.12 0.03 0.86 -selax \
	-save figs/eos_s2.png -show &

splot-eos:
	o2graph -subplots 4 1 sharex=True \
	-subadj hspace=0.0,left=0.16,right=0.8,top=0.99,bottom=0.05 \
	-selax 0 -set logy 1 -read files/esmr_NL.o2 ed_pr_3d \
	-den-plot prob pcm=True,cmap=pink_r \
	-read files/esmr_NL.o2 c68_es -plot lw=1,color=purple \
	-read files/esmr_NL.o2 c95_es -plot lw=1,color=orange \
	-ttext 0.1 0.92 "NL" \
	-selax 1 -set logy 1 -read files/esmr_NP.o2 ed_pr_3d \
	-den-plot prob pcm=True,cmap=pink_r \
	-read files/esmr_NP.o2 c68_es -plot lw=1,color=purple \
	-read files/esmr_NP.o2 c95_es -plot lw=1,color=orange \
	-ttext 0.1 0.92 "NP" \
	-selax 2 -set logy 1 -read files/esmr_ML.o2 ed_pr_3d \
	-den-plot prob pcm=True,cmap=pink_r \
	-read files/esmr_ML.o2 c68_es -plot lw=1,color=purple \
	-read files/esmr_ML.o2 c95_es -plot lw=1,color=orange \
	-ttext 0.1 0.92 "ML" \
	-ytitle ".$$ \qquad \qquad \qquad \qquad \qquad \qquad \quad $$ Pressure, $$ P~[\mathrm{MeV}/\mathrm{fm}^3] $$" \
	-selax 3 -set logy 1 -read files/esmr_MP.o2 ed_pr_3d \
	-den-plot prob pcm=True,cmap=pink_r \
	-read files/esmr_MP.o2 c68_es -plot lw=1,color=purple \
	-read files/esmr_MP.o2 c95_es -plot lw=1,color=orange \
	-ttext 0.1 0.92 "MP" \
	-xtitle "Energy Density, $$ \epsilon~[\mathrm{MeV}/\mathrm{fm}^3] $$" \
	-addcbar 0.82 0.05 0.03 0.94 -selax \
	-show &

splot-mvsr1:
	o2graph -subplots 2 1 sharex=True \
	-subadj hspace=0.0,left=0.14,right=0.84,top=0.98,bottom=0.12 \
	-selax 0 -read files/esmr_NL.o2 gm_r_3d -den-plot prob cmap=pink_r \
	-read files/esmr_NL.o2 c68_mr -plot lw=1,color=purple \
	-read files/esmr_NL.o2 c95_mr -plot lw=1,color=orange \
	-ttext 0.1 0.92 "NL" \
	-selax 1 -read files/esmr_NP.o2 gm_r_3d -den-plot prob cmap=pink_r \
	-read files/esmr_NP.o2 c68_mr -plot lw=1,color=purple \
	-read files/esmr_NP.o2 c95_mr -plot lw=1,color=orange \
	-ttext 0.1 0.92 "NP" \
	-xtitle "Radius, $$ R~[\mathrm{km}] $$" \
	-ytitle ".$$ \qquad \qquad \qquad \qquad \qquad \quad $$ Mass, $$ M~[\mathrm{M}_\odot] $$" \
	-addcbar 0.86 0.12 0.03 0.86 -selax \
	-save figs/mvsr_s1.png -show &

splot-mvsr2:
	o2graph -subplots 2 1 sharex=True \
	-subadj hspace=0.0,left=0.14,right=0.84,top=0.98,bottom=0.12 \
	-selax 0 -read files/esmr_ML.o2 gm_r_3d -den-plot prob cmap=pink_r \
	-read files/esmr_ML.o2 c68_mr -plot lw=1,color=purple \
	-read files/esmr_ML.o2 c95_mr -plot lw=1,color=orange \
	-ttext 0.1 0.92 "ML" \
	-selax 1 -read files/esmr_MP.o2 gm_r_3d -den-plot prob cmap=pink_r \
	-read files/esmr_MP.o2 c68_mr -plot lw=1,color=purple \
	-read files/esmr_MP.o2 c95_mr -plot lw=1,color=orange \
	-xtitle "Radius, $$ R~[\mathrm{km}] $$" \
	-ytitle ".$$ \qquad \qquad \qquad \qquad \qquad \quad $$ Mass, $$ M~[\mathrm{M}_\odot] $$" \
	-ttext 0.1 0.92 "MP" \
	-addcbar 0.86 0.12 0.03 0.86 -selax \
	-save figs/mvsr_s2.png -show &

splot-mvsr:
	o2graph -subplots 4 1 sharex=True \
	-subadj hspace=0.0,left=0.16,right=0.8,top=0.99,bottom=0.05 \
	-selax 0 -read files/esmr_NL.o2 gm_r_3d \
	-den-plot prob cmap=pink_r \
	-read files/esmr_NL.o2 c68_mr -plot lw=1,color=purple \
	-read files/esmr_NL.o2 c95_mr -plot lw=1,color=orange \
	-ttext 0.1 0.92 "NL" \
	-selax 1 -read files/esmr_NP.o2 gm_r_3d \
	-den-plot prob cmap=pink_r \
	-read files/esmr_NP.o2 c68_mr -plot lw=1,color=purple \
	-read files/esmr_NP.o2 c95_mr -plot lw=1,color=orange \
	-ttext 0.1 0.92 "NP" \
	-selax 2 -read files/esmr_ML.o2 gm_r_3d \
	-den-plot prob cmap=pink_r \
	-read files/esmr_ML.o2 c68_mr -plot lw=1,color=purple \
	-read files/esmr_ML.o2 c95_mr -plot lw=1,color=orange \
	-ttext 0.1 0.92 "ML" \
	-ytitle ".$$ \qquad \qquad \qquad \qquad \qquad \qquad \quad $$ Mass, $$ M~[\mathrm{M}_\odot] $$" \
	-selax 3 -read files/esmr_MP.o2 gm_r_3d \
	-den-plot prob cmap=pink_r \
	-read files/esmr_MP.o2 c68_mr -plot lw=1,color=purple \
	-read files/esmr_MP.o2 c95_mr -plot lw=1,color=orange \
	-ttext 0.1 0.92 "MP" \
	-xtitle "Radius, $$ R~[\mathrm{km}] $$" \
	-addcbar 0.82 0.05 0.03 0.94 -selax \
	-show &

plot-hist-mmax:
	o2graph -read files/main2_$(model).o2 markov_chain_0 \
	-to-hist M_max 40 log_wgt -plot density=True \
	-xtitle "M_max [$$ M_\odot $$]" -ytitle "Frequency" \
	-ttext 0.5 0.96 "$(model): $(title)" \
	-save figs/hist_mmax_$(model).png -show &

plot-scat-r-sn:
	o2graph -read files/main2_$(model).o2 markov_chain_0 -function 1 w -set colbar 1 \
	-scatter R_43 SN_NS w log_wgt cmap=cividis \
	-xtitle "Radius $$ R_{1.4M_\odot} $$ [km]" \
	-ytitle "Mass Distribution, $$ \mathcal{P}(M=1.4 M_\odot $$)" \
	-ttext 0.5 0.96 "$(model): DNS" -save figs/scat_sn-r_ns_$(model).png -show &
	o2graph -read files/main2_$(model).o2 markov_chain_0 -function 1 w -set colbar 1 \
	-scatter R_50 SN_WD w log_wgt cmap=cividis \
	-xtitle "Radius $$ R_{1.6M_\odot} $$ [km]" \
	-ytitle "Mass Distribution, $$ \mathcal{P}(M=1.6 M_\odot $$)" \
	-ttext 0.5 0.96 "$(model): NS-WD" -save figs/scat_sn-r_wd_$(model).png -show &
	o2graph -read files/main2_$(model).o2 markov_chain_0 -function 1 w -set colbar 1 \
	-scatter R_53 SN_LM w log_wgt cmap=cividis \
	-xtitle "Radius $$ R_{1.7M_\odot} $$ [km]" \
	-ytitle "Mass Distribution, $$ \mathcal{P}(M=1.7 M_\odot $$)" \
	-ttext 0.5 0.96 "$(model): LMXB" -save figs/scat_sn-r_lm_$(model).png -show &

splot-r-sn:
	o2graph -read files/main2_$(model).o2 markov_chain_0 -subplots 3 1 \
	-subadj hspace=0.8,left=0.12,right=0.8,top=0.94,bottom=0.12 \
	-selax 0 -function 1 w -scatter R_43 SN_NS w log_wgt cmap=cividis \
	-xtitle "$$ R_{1.4M_\odot} $$ [km]" \
	-ytitle "$$ \mathcal{P}(M=1.4 M_\odot $$)" \
	-ttext 0.5 1.1 "$(model): DNS" \
	-selax 1 -function 1 w -scatter R_50 SN_WD w log_wgt cmap=cividis \
	-xtitle "$$ R_{1.6M_\odot} $$ [km]" \
	-ytitle "$$ \mathcal{P}(M=1.6 M_\odot $$)" \
	-ttext 0.5 1.1 "$(model): NS-WD" \
	-selax 2 -function 1 w -scatter R_53 SN_LM w log_wgt cmap=cividis \
	-xtitle "$$ R_{1.7M_\odot} $$ [km]" \
	-ytitle "$$ \mathcal{P}(M=1.7 M_\odot $$)" \
	-ttext 0.5 1.1 "$(model): LMXB" \
	-addcbar 0.82 0.12 0.04 0.82 -selax \
	-save figs/sn_r_$(model).png -show &

splot-kde-lmxb:
	o2graph -read files/main2_NL.o2 markov_chain_0 \
	-kde-2d-plot R_53 SN_LM \
	-ttext 0.1 0.92 "NL" -save figs/scat-r53_NL.pdf -show &
	o2graph -read files/main2_NL.o2 markov_chain_0 \
	-kde-2d-plot Phc_23 SN_LM \
	-save figs/scat-p23_NL.pdf -show &
	o2graph -read files/main2_NL.o2 markov_chain_0 \
	-kde-2d-plot M_max SN_LM \
	-save figs/scat-mmax_NL.pdf-show &
	o2graph -read files/main2_NP.o2 markov_chain_0 \
	-kde-2d-plot R_53 SN_LM \
	-ttext 0.1 0.92 "NP" \
	-ytitle "$$ \mathcal{P}(M=1.7 M_\odot $$)" -show &
	o2graph -read files/main2_NP.o2 markov_chain_0 \
	-kde-2d-plot Phc_23 SN_LM -show &
	o2graph -read files/main2_NP.o2 markov_chain_0 \
	-kde-2d-plot M_max SN_LM -show &	
	o2graph -read files/main2_ML.o2 markov_chain_0 \
	-kde-2d-plot R_53 SN_LM \
	-ttext 0.1 0.92 "ML" -show &
	o2graph -read files/main2_ML.o2 markov_chain_0 \
	-kde-2d-plot Phc_23 SN_LM -show &
	o2graph -read files/main2_ML.o2 markov_chain_0 \
	-kde-2d-plot M_max SN_LM -show &
	o2graph -read files/main2_MP.o2 markov_chain_0 \
	-kde-2d-plot R_53 SN_LM \
	-ttext 0.1 0.92 "MP" \
	-xtitle "$$ R_{1.7M_\odot} $$ [km]" -show &
	o2graph -read files/main2_MP.o2 markov_chain_0 \
	-kde-2d-plot Phc_23 SN_LM \
	-xtitle "$$ P~(\epsilon \approx $$ 500 MeV/fm $$ ^3) $$" -show &
	o2graph -read files/main2_MP.o2 markov_chain_0 \
	-kde-2d-plot M_max SN_LM \
	-xtitle "$$ M_\mathrm{max} [M_\odot] $$" -show &

plot-scat-p23-sn:
	o2graph -read files/main2_$(model).o2 markov_chain_0 -function 1 w -set colbar 1 \
	-scatter Phc_23 SN_NS w log_wgt cmap=cividis \
	-xtitle "Pressure $$ P~(\epsilon \approx $$ 500 MeV/fm $$ ^3) $$" \
	-ytitle "Mass Distribution, $$ \mathcal{P}(M=1.4 M_\odot $$)" \
	-ttext 0.5 0.96 "$(model): DNS" \
	-save figs/scat_sn-p23_$(model).png -show &
	o2graph -read files/main2_$(model).o2 markov_chain_0 -function 1 w -set colbar 1 \
	-scatter Phc_23 SN_WD w log_wgt cmap=cividis \
	-xtitle "Pressure $$ P~(\epsilon \approx $$ 500 MeV/fm $$ ^3) $$" \
	-ytitle "Mass Distribution, $$ \mathcal{P}(M=1.6 M_\odot $$)" \
	-ttext 0.5 0.96 "$(model): NS-WD" \
	-save figs/scat_sn-p23_$(model).png -show &
	o2graph -read files/main2_$(model).o2 markov_chain_0 -function 1 w -set colbar 1 \
	-scatter Phc_23 SN_LM w log_wgt cmap=cividis \
	-xtitle "Pressure $$ P~(\epsilon \approx $$ 500 MeV/fm $$ ^3) $$" \
	-ytitle "Mass Distribution, $$ \mathcal{P}(M=1.7 M_\odot $$)" \
	-ttext 0.5 0.96 "$(model): LMXB" \
	-save figs/scat_sn-p23_$(model).png -show &

splot-p23-sn:
	o2graph -read files/main2_$(model).o2 markov_chain_0 -subplots 3 1 sharex=True \
	-subadj hspace=0.6,left=0.12,right=0.8,top=0.94,bottom=0.12 \
	-selax 0 -function 1 w -scatter Phc_23 SN_NS w log_wgt cmap=cividis \
	-ttext 0.5 1.1 "$(model): DNS" \
	-ytitle "$$ \mathcal{P}(M=1.4 M_\odot $$)" \
	-selax 1 -function 1 w -scatter Phc_23 SN_WD w log_wgt cmap=cividis \
	-ttext 0.5 1.1 "$(model): NS-WD" \
	-ytitle "$$ \mathcal{P}(M=1.6 M_\odot $$)" \
	-selax 2 -function 1 w -scatter Phc_23 SN_LM w log_wgt cmap=cividis \
	-ttext 0.5 1.1 "$(model): LMXB" \
	-xtitle "$$ P~(\epsilon \approx $$ 500 MeV/fm $$ ^3) $$" \
	-ytitle "$$ \mathcal{P}(M=1.7 M_\odot $$)" \
	-addcbar 0.82 0.12 0.04 0.82 -selax \
	-save figs/scat_sn-p23_$(model).png -show &

plot-scat-mmax-sn:
	o2graph -read files/main2_$(model).o2 markov_chain_0 -function 1 w -set colbar 1 \
	-scatter M_max SN_NS w log_wgt cmap=cividis \
	-xtitle "Maximum Mass $$ M_\mathrm{max} [M_\odot] $$" \
	-ytitle "Mass Distribution, $$ \mathcal{P}(M=1.4 M_\odot $$)" \
	-ttext 0.5 0.96 "$(model): DNS" \
	-save figs/scat_sn-mmax_$(model).png -show &
	o2graph -read files/main2_$(model).o2 markov_chain_0 -function 1 w -set colbar 1 \
	-scatter M_max SN_WD w log_wgt cmap=cividis \
	-xtitle "Maximum Mass $$ M_\mathrm{max} [M_\odot] $$" \
	-ytitle "Mass Distribution, $$ \mathcal{P}(M=1.6 M_\odot $$)" \
	-ttext 0.5 0.96 "$(model): NS-WD" \
	-save figs/scat_sn-mmax_$(model).png -show &
	o2graph -read files/main2_$(model).o2 markov_chain_0 -function 1 w -set colbar 1 \
	-scatter M_max SN_LM w log_wgt cmap=cividis \
	-xtitle "Maximum Mass $$ M_\mathrm{max} [M_\odot] $$" \
	-ytitle "Mass Distribution, $$ \mathcal{P}(M=1.7 M_\odot $$)" \
	-ttext 0.5 0.96 "$(model): LMXB" \
	-save figs/scat_sn-mmax_$(model).png -show &

splot-mmax-sn:
	o2graph -read files/main2_$(model).o2 markov_chain_0 -subplots 3 1 sharex=True \
	-subadj hspace=0.6,left=0.12,right=0.8,top=0.94,bottom=0.12 \
	-selax 0 -function 1 w -scatter M_max SN_NS w log_wgt cmap=cividis \
	-ttext 0.5 1.1 "$(model): DNS" \
	-ytitle "$$ \mathcal{P}(M=1.4 M_\odot $$)" \
	-selax 1 -function 1 w -scatter M_max SN_WD w log_wgt cmap=cividis \
	-ttext 0.5 1.1 "$(model): NS-WD" \
	-ytitle "$$ \mathcal{P}(M=1.6 M_\odot $$)" \
	-selax 2 -function 1 w -scatter M_max SN_LM w log_wgt cmap=cividis \
	-ttext 0.5 1.1 "$(model): LMXB" \
	-xtitle "$$ M_\mathrm{max} [M_\odot] $$" \
	-ytitle "$$ \mathcal{P}(M=1.7 M_\odot $$)" \
	-addcbar 0.82 0.12 0.04 0.82 -selax \
	-save figs/scat_sn-p23_$(model).png -show &

splot-dist1:
	o2graph -subplots 2 3 sharex=True,sharey=True \
	-subadj hspace=0.0,left=0.05,right=0.93,top=0.95,bottom=0.12 \
	-selax 0 -read files/dist_NL.o2 sn_ns -den-plot prob_ns cmap=pink_r \
	-read files/dist_NL.o2 c68_ns -plot lw=1,color=purple \
	-read files/dist_NL.o2 c95_ns -plot lw=1,color=orange \
	-ttext 0.5 1.05 "DNS" \
	-ttext 0.2 0.92 "NL" \
	-ytitle "Probability Density, $$ \mathcal{P}(M) $$" \
	-selax 1 -read files/dist_NL.o2 sn_wd -den-plot prob_wd cmap=pink_r \
	-read files/dist_NL.o2 c68_wd -plot lw=1,color=purple \
	-read files/dist_NL.o2 c95_wd -plot lw=1,color=orange \
	-ttext 0.5 1.05 "NS-WD" \
	-selax 2 -read files/dist_NL.o2 sn_lm -den-plot prob_lm cmap=pink_r \
	-read files/dist_NL.o2 c68_lm -plot lw=1,color=purple \
	-read files/dist_NL.o2 c95_lm -plot lw=1,color=orange \
	-ttext 0.5 1.05 "LMXB" \
	-selax 3 -read files/dist_NP.o2 sn_ns -den-plot prob_ns cmap=pink_r \
	-read files/dist_NP.o2 c68_ns -plot lw=1,color=purple \
	-read files/dist_NP.o2 c95_ns -plot lw=1,color=orange \
	-ttext 0.2 0.92 "NP" \
	-selax 4 -read files/dist_NP.o2 sn_wd -den-plot prob_wd cmap=pink_r \
	-read files/dist_NP.o2 c68_wd -plot lw=1,color=purple \
	-read files/dist_NP.o2 c95_wd -plot lw=1,color=orange \
	-xtitle "Mass, $$ M $$ [$$ M_\odot $$]" \
	-selax 5 -read files/dist_NP.o2 sn_lm -den-plot prob_lm cmap=pink_r \
	-read files/dist_NP.o2 c68_lm -plot lw=1,color=purple \
	-read files/dist_NP.o2 c95_lm -plot lw=1,color=orange \
	-addcbar 0.94 0.12 0.02 0.83 -selax \
	-save figs/sn_s1.png -show &

splot-dist2:
	o2graph -subplots 2 3 sharex=True,sharey=True \
	-subadj hspace=0.0,left=0.05,right=0.93,top=0.95,bottom=0.12 \
	-selax 0 -read files/dist_ML.o2 sn_ns -den-plot prob_ns cmap=pink_r \
	-read files/dist_ML.o2 c68_ns -plot lw=1,color=purple \
	-read files/dist_ML.o2 c95_ns -plot lw=1,color=orange \
	-ttext 0.5 1.05 "DNS" \
	-ttext 0.2 0.92 "ML" \
	-ytitle "Probability Density, $$ \mathcal{P}(M) $$" \
	-selax 1 -read files/dist_ML.o2 sn_wd -den-plot prob_wd cmap=pink_r \
	-read files/dist_ML.o2 c68_wd -plot lw=1,color=purple \
	-read files/dist_ML.o2 c95_wd -plot lw=1,color=orange \
	-ttext 0.5 1.05 "NS-WD" \
	-selax 2 -read files/dist_ML.o2 sn_lm -den-plot prob_lm cmap=pink_r \
	-read files/dist_ML.o2 c68_lm -plot lw=1,color=purple \
	-read files/dist_ML.o2 c95_lm -plot lw=1,color=orange \
	-ttext 0.5 1.05 "LMXB" \
	-selax 3 -read files/dist_MP.o2 sn_ns -den-plot prob_ns cmap=pink_r \
	-read files/dist_MP.o2 c68_ns -plot lw=1,color=purple \
	-read files/dist_MP.o2 c95_ns -plot lw=1,color=orange \
	-ttext 0.2 0.92 "MP" \
	-selax 4 -read files/dist_MP.o2 sn_wd -den-plot prob_wd cmap=pink_r \
	-read files/dist_MP.o2 c68_wd -plot lw=1,color=purple \
	-read files/dist_MP.o2 c95_wd -plot lw=1,color=orange \
	-xtitle "Mass, $$ M $$ [$$ M_\odot $$]" \
	-selax 5 -read files/dist_MP.o2 sn_lm -den-plot prob_lm cmap=pink_r \
	-read files/dist_MP.o2 c68_lm -plot lw=1,color=purple \
	-read files/dist_MP.o2 c95_lm -plot lw=1,color=orange \
	-addcbar 0.94 0.12 0.02 0.83 -selax \
	-save figs/sn_s2.png -show &

splot-dist:
	o2graph -subplots 4 3 sharex=True,sharey=True \
	-subadj hspace=0.0,left=0.06,right=0.92,top=0.98,bottom=0.05 \
	-selax 0 -read files/dist_NL.o2 sn_ns -den-plot prob_ns cmap=pink_r \
	-read files/dist_NL.o2 c68_ns -plot lw=1,color=purple \
	-read files/dist_NL.o2 c95_ns -plot lw=1,color=orange \
	-ttext 0.2 0.92 "NL" \
	-ttext 0.5 1.05 "DNS" \
	-selax 1 -read files/dist_NL.o2 sn_wd -den-plot prob_wd cmap=pink_r \
	-read files/dist_NL.o2 c68_wd -plot lw=1,color=purple \
	-read files/dist_NL.o2 c95_wd -plot lw=1,color=orange \
	-ttext 0.5 1.05 "NS-WD" \
	-selax 2 -read files/dist_NL.o2 sn_lm -den-plot prob_lm cmap=pink_r \
	-read files/dist_NL.o2 c68_lm -plot lw=1,color=purple \
	-read files/dist_NL.o2 c95_lm -plot lw=1,color=orange \
	-ttext 0.5 1.05 "LMXB" \
	-selax 3 -read files/dist_NP.o2 sn_ns -den-plot prob_ns cmap=pink_r \
	-read files/dist_NP.o2 c68_ns -plot lw=1,color=purple \
	-read files/dist_NP.o2 c95_ns -plot lw=1,color=orange \
	-ttext 0.2 0.92 "NP" \
	-selax 4 -read files/dist_NP.o2 sn_wd -den-plot prob_wd cmap=pink_r \
	-read files/dist_NP.o2 c68_wd -plot lw=1,color=purple \
	-read files/dist_NP.o2 c95_wd -plot lw=1,color=orange \
	-xtitle "Mass, $$ M $$ [$$ M_\odot $$]" \
	-selax 5 -read files/dist_NP.o2 sn_lm -den-plot prob_lm cmap=pink_r \
	-read files/dist_NP.o2 c68_lm -plot lw=1,color=purple \
	-read files/dist_NP.o2 c95_lm -plot lw=1,color=orange \
	-selax 6 -read files/dist_ML.o2 sn_ns -den-plot prob_ns cmap=pink_r \
	-read files/dist_ML.o2 c68_ns -plot lw=1,color=purple \
	-read files/dist_ML.o2 c95_ns -plot lw=1,color=orange \
	-ttext 0.2 0.92 "ML" \
	-ytitle ".$$ \qquad \qquad \qquad \qquad \qquad \qquad \quad $$ Probability Density, $$ \mathcal{P}(M) $$" \
	-selax 7 -read files/dist_ML.o2 sn_wd -den-plot prob_wd cmap=pink_r \
	-read files/dist_ML.o2 c68_wd -plot lw=1,color=purple \
	-read files/dist_ML.o2 c95_wd -plot lw=1,color=orange \
	-selax 8 -read files/dist_ML.o2 sn_lm -den-plot prob_lm cmap=pink_r \
	-read files/dist_ML.o2 c68_lm -plot lw=1,color=purple \
	-read files/dist_ML.o2 c95_lm -plot lw=1,color=orange \
	-selax 9 -read files/dist_MP.o2 sn_ns -den-plot prob_ns cmap=pink_r \
	-read files/dist_MP.o2 c68_ns -plot lw=1,color=purple \
	-read files/dist_MP.o2 c95_ns -plot lw=1,color=orange \
	-ttext 0.2 0.92 "MP" \
	-selax 10 -read files/dist_MP.o2 sn_wd -den-plot prob_wd cmap=pink_r \
	-read files/dist_MP.o2 c68_wd -plot lw=1,color=purple \
	-read files/dist_MP.o2 c95_wd -plot lw=1,color=orange \
	-xtitle "Mass, $$ M $$ [$$ M_\odot $$]" \
	-selax 11 -read files/dist_MP.o2 sn_lm -den-plot prob_lm cmap=pink_r \
	-read files/dist_MP.o2 c68_lm -plot lw=1,color=purple \
	-read files/dist_MP.o2 c95_lm -plot lw=1,color=orange \
	-addcbar 0.93 0.05 0.02 0.93 -selax \
	-save figs/sn_s.png -show &

file-eos:
	acol -read files/main_$(model).o2 markov_chain_0 \
	$(foreach i,$(INDEX),-function P_$i*197.32698 Phc_$i) \
	-set obj_name ed_pr -select ":Phc_*" -internal files/eos_$(model).o2 \
	-read files/eos_$(model).o2 ed_pr -ser-hist-t3d \
	"grid:0.3*197.32698,10.0*197.32698,((10.0-0.3)*197.32698)/99.0" x ed \
	"grid:1.0,2.0e3,(2.0e3)^(1.0/99.0),log" auto pr \
	"Phc_*" prob -type -list \
	-set obj_name ed_pr_3d -internal files/eos_$(model).o2 \
	-to-hist-2d prob \
	-contours frac 0.68 files/eos_$(model).o2 c68 \
	-contours frac 0.95 files/eos_$(model).o2 c95

plot-eos:
	o2graph -read files/eos_$(model).o2 ed_pr_3d -set logy 1 -set colbar 1 \
	-den-plot prob pcm=True,cmap=pink_r \
	-read files/eos_$(model).o2 c68 -plot lw=1,color=purple \
	-read files/eos_$(model).o2 c95 -plot lw=1,color=orange \
	-xtitle "Energy Density, $$ \epsilon~[\mathrm{MeV}/\mathrm{fm}^3] $$" \
	-ytitle "Pressure, $$ P~[\mathrm{MeV}/\mathrm{fm}^3] $$" \
	-ttext 0.5 0.96 "$(model): EoS" \
	-save figs/eos_$(model).png -show &

file-mvsr:
	acol -read files/main_$(model).o2 markov_chain_0 -delete-col R_max \
	-set obj_name gm_r -select ":R_*" -internal files/mvsr_$(model).o2 \
	-read files/mvsr_$(model).o2 gm_r -ser-hist-t3d \
	"grid:0.2,3.0,(3.0-0.2)/99.0" y gm \
	"grid:9.0,15.0,(15.0-9.0)/99.0" auto r \
	"R_*" prob -type -list \
	-set obj_name gm_r_3d -internal files/mvsr_$(model).o2 \
	-to-hist-2d prob \
	-contours frac 0.68 files/mvsr_$(model).o2 c68 \
	-contours frac 0.95 files/mvsr_$(model).o2 c95

plot-mvsr:
	o2graph -read files/mvsr_$(model).o2 gm_r_3d -set colbar 1 \
	-den-plot prob cmap=pink_r \
	-read files/mvsr_$(model).o2 c68 -plot lw=1,color=purple \
	-read files/mvsr_$(model).o2 c95 -plot lw=1,color=orange \
	-xtitle "Radius, $$ R~[\mathrm{km}] $$" \
	-ytitle "Mass, $$ M~[\mathrm{M}_\odot] $$" \
	-ttext 0.5 0.96 "$(model): $$ M-R $$" \
	-save figs/mvsr_$(model).png -exec plots.py -show &

file-eoa:
	acol -read files/main_$(model).o2 markov_chain_0 -select ":EoA_*" \
	-set obj_name nb_ea -internal files/eoa_$(model).o2 \
	-read files/eoa_$(model).o2 nb_ea -ser-hist-t3d \
	"grid:0.04,1.24,(1.24-0.04)/99.0" x nb \
	"grid:0.0,2.5,(2.5-0.0)/99.0" auto ea \
	"EoA_*" prob -type -list \
	-set obj_name nb_ea_3d -internal files/eoa_$(model).o2 \
	-to-hist-2d prob \
	-contours frac 0.68 files/eoa_$(model).o2 c68 \
	-contours frac 0.95 files/eoa_$(model).o2 c95 \
	-contours frac 0.99 files/eoa_$(model).o2 c99

plot-eoa:
	o2graph -read files/eoa_$(model).o2 nb_ea_3d -set colbar 1 \
	-den-plot prob cmap=pink_r \
	-read files/eoa_$(model).o2 c68 -plot lw=1,color=purple \
	-read files/eoa_$(model).o2 c95 -plot lw=1,color=orange \
	-xtitle "Baryon Density, $$ n_b~[\mathrm{fm}^{-3}] $$" \
	-ytitle "Energy/Nucleon, $$ E/A~[\mathrm{MeV}] $$" \
	-ttext 0.5 0.96 "$(model): $(title)" \
	-save figs/eoa_$(model).png -show &

file-pnb:
	acol -read files/main_$(model).o2 markov_chain_0 \
	$(foreach i,$(INDEX),-function Pnb_$i*197.32698 Pnb_hc_$i) \
	-select ":Pnb_hc_*" -set obj_name nb_pr -internal files/pnb_$(model).o2 \
	-read files/pnb_$(model).o2 nb_pr -ser-hist-t3d \
	"grid:0.04,1.24,(1.24-0.04)/99.0" x nb \
	"grid:1.0,2.0e3,(2.0e3)^(1/99),log" auto pr \
	"Pnb_hc_*" prob -type -list \
	-set obj_name nb_pr_3d -internal files/pnb_$(model).o2 \
	-to-hist-2d prob \
	-contours frac 0.68 files/pnb_$(model).o2 c68 \
	-contours frac 0.95 files/pnb_$(model).o2 c95 \
	-contours frac 0.99 files/pnb_$(model).o2 c99

plot-pnb:
	o2graph -read files/pnb_$(model).o2 nb_pr_3d \
	-set colbar 1 -set logy 1 \
	-den-plot prob pcm=True,cmap=pink_r \
	-read files/pnb_$(model).o2 c68 -plot lw=1,color=purple \
	-read files/pnb_$(model).o2 c95 -plot lw=1,color=orange \
	-xtitle "Baryon Density, $$ n_b~[\mathrm{fm}^{-3}] $$" \
	-ytitle "Pressure, $$ P~[\mathrm{MeV/fm}^3] $$" \
	-ttext 0.5 0.96 "$(model): $(title)" \
	-save figs/pnb_$(model).png -show &

plot-pnb2:
	o2graph -plotv "grid:0.04,1.24,(1.24-0.04)/99.0" \
	"hdf5:files/pnb_ML.o2:nb_pr:0-99:Pnb_hc_*" \
	-xtitle "Baryon Density, $ n_b~[\mathrm{fm}^{-3}] $" \
	-ytitle "Pressure, $ P~[\mathrm{MeV/fm}^3] $" -show &

plot-trc-wgt:
	o2graph -read files/main_$(model).o2 markov_chain_0 \
	-xtitle "Steps" -ytitle "Average log_wgt" \
	-ttext 0.5 0.96 "$(model): $(title)" \
	-average-rows log_wgt 512 -plot1 log_wgt \
	-save figs/trc_wgt_$(model).png -show &

plot-trc-mmax:
	o2graph -read files/main_$(model).o2 markov_chain_0 \
	-xtitle "Steps" -ytitle "Average M_max" \
	-ttext 0.5 0.96 "$(model): $(title)" \
	-average-rows M_max 512 -plot1 M_max \
	-save figs/trc_mmax_$(model).png -show &

plot-trc-a:
	o2graph -read files/main_$(model).o2 markov_chain_0 \
	-xtitle "Steps" -ytitle "Average $$ a $$ [MeV]" \
	-ttext 0.5 0.96 "$(model): $(title)" \
	-average-rows a 512 -plot1 a \
	-save figs/trc_a_$(model).png -show &

plot-trc-alpha:
	o2graph -read files/main_$(model).o2 markov_chain_0 \
	-xtitle "Steps" -ytitle "Average $$ \alpha $$ [MeV]" \
	-ttext 0.5 0.96 "$(model): $(title)" \
	-average-rows alpha 512 -plot1 alpha \
	-save figs/trc_alpha_$(model).png -show &

plot-trc-r43:
	o2graph -read files/main_$(model).o2 markov_chain_0 \
	-xtitle "Steps" -ytitle "Average $$ R_{1.4M_\odot} $$ [km]" \
	-ttext 0.5 0.96 "$(model): $(title)" \
	-average-rows R_43 512 -plot1 R_43 \
	-save figs/trc_r43_$(model).png -show &

plot-trc-dpdm:
	o2graph -read files/main_$(model).o2 markov_chain_0 \
	-function "log(dpdM)" log_dpdm \
	-xtitle "Steps" -ytitle "Average log(dp/dM)" \
	-ttext 0.5 0.96 "$(model): $(title)" \
	-average-rows log_dpdm 512 -plot1 log_dpdm \
	-save figs/trc_dpdm_$(model).png -show &

plot-acor-a:
	o2graph $(foreach i, $(THREAD),$(foreach j,$(WALKER), \
	-read files/main_$(model).o2 markov_chain_0 \
	-select-rows "thread==$i && walker==$j" \
	-autocorr store,default a -plot 0)) \
	-xtitle "Autocorr. Length" -ytitle "Autocorr. Coeff." \
	-ttext 0.5 0.96 "$(model): $(title)" \
	-save figs/a-acor_$(model).png -show &
	
file-all:
	make model="ML" file-esmr  && make model="MP" file-esmr  && make model="NL" file-esmr  && make model="NP" file-esmr  && \
	make model="ML" file-main2 && make model="MP" file-main2 && make model="NL" file-main2 && make model="NP" file-main2 && \
	make model="ML" file-dist  && make model="MP" file-dist  && make model="NL" file-dist  && make model="NP" file-dist

fname =

file-min:
	acol -read $(fname) markov_chain_0 -select rank thread walker mult log_wgt a alpha param_S param_L csq1 trans1 csq2 trans2 csq3 M_chirp_det q z_cdf m1_gw19 mf_6304 mf_6397 mf_M13 mf_M28 mf_M30 mf_wCen mf_X7 mf_1810b mf_1724b mf_1702 mf_0030 mf_0740 mean_NS log10_width_NS skewness_NS mean_WD log10_width_WD skewness_WD mean_LMS log10_width_LMS skewness_LMS M_J0453p M_J0453c M_J1906p M_J1906c M_B1534p M_B1534c M_B1913p M_B1913c M_B2127p M_B2127c M_J0737A M_J0737B M_J1756p M_J1756c M_J1807p M_J1807c M_J1518p M_J1518c M_J1811p M_J1811c M_J1829p M_J1829c M_J2045 M_J2053 M_J1713 M_B1855 M_J0751 M_J1141 M_J1738 M_J1614 M_J0348 M_J2222 M_J2234 M_J1949 M_J1012 M_J0437 M_J1909 M_J1802 M_J1911 M_J2043 M_J0337 M_J1946 M_J1918 M_J1600 M_J0621 M_B2303 M_J0024 M_J0514 M_B1516 M_J1748I M_J1748J M_B1802 M_B1911 M_J0740 M_CygX2 M_XTEJ2123 M_4U1822 M_HerX1 M_2S0921 -internal $(fname)_in

models = ml mp nl np
rank = $(shell seq 0 15)
rank2 = $(shell seq 1 15)
rank3 = 0 1 2 3 4 5 6 7 8 9 11 12 13 14 15

scp-down-br2:
	scp manik@br2:/ocean/projects/phy230011p/shared/bamr/out/*_out .

file-mult:
	$(foreach m,$(models),\
		$(foreach r,$(rank),\
			acol -read $(m)_$(r)_out -select-rows "mult>0.5" -internal $(m)_$(r)_out;\
		)\
	)

file-hcp:
	$(foreach m,$(models),\
		$(foreach r,$(rank),\
			$(if $(filter $(m),ml),\
				acol -h5-copy $(m)_$(r)_out $(m)_65a$(r);,\
				$(if $(filter $(m),mp),\
					acol -h5-copy $(m)_$(r)_out $(m)_56a$(r);,\
					acol -h5-copy $(m)_$(r)_out $(m)_54a$(r);\
				)\
			)\
		)\
	)

file-con:
	$(foreach m,$(models),\
		$(if $(filter $(m),ml),\
			acol -read $(m)_65a0 $(foreach r,$(rank2),-cat $(m)_65a$(r)) -internal $(m)_65;,\
			$(if $(filter $(m),mp),\
				acol -read $(m)_56a0 $(foreach r,$(rank2),-cat $(m)_56a$(r)) -internal $(m)_56;,\
				acol -read $(m)_54a0 $(foreach r,$(rank2),-cat $(m)_54a$(r)) -internal $(m)_54;\
			)\
		)\
	)

scp-up-br2:
	scp ml_65a* manik@br2:/ocean/projects/phy230011p/shared/bamr/out/.
	scp mp_56a* manik@br2:/ocean/projects/phy230011p/shared/bamr/out/.
	scp nl_54a* manik@br2:/ocean/projects/phy230011p/shared/bamr/out/.
	scp np_54a* manik@br2:/ocean/projects/phy230011p/shared/bamr/out/.

scp-up-iso:
	scp ml_65a* anikmh@iso:/mnt/hyperon/anikmh/bridges/r65_16n_8c_72h/.
	scp mp_56a* anikmh@iso:/mnt/hyperon/anikmh/bridges/r56_16n_8c_72h/.
	scp nl_54a* anikmh@iso:/mnt/hyperon/anikmh/bridges/r54_16n_8c_72h/.
	scp np_54a* anikmh@iso:/mnt/hyperon/anikmh/bridges/r54_16n_8c_72h/.

ML = $(shell seq 32 63)
MP = $(shell seq 33 54)
NL = $(shell seq 29 51)
NP = $(shell seq 31 51)

file-nlines:
	acol -read ml_31 markov_chain_0 $(foreach i,$(ML),-cat ml_$i) -nlines
	acol -read mp_32 markov_chain_0 $(foreach i,$(MP),-cat mp_$i) -nlines
	acol -read nl_28 markov_chain_0 $(foreach i,$(NL),-cat nl_$i) -nlines
	acol -read np_30 markov_chain_0 $(foreach i,$(NP),-cat np_$i) -nlines

file-rm-last:
	$(foreach i,$(ML),acol -read ml_$i markov_chain_0 -index -select-rows "N>511" -delete-col N -internal files/ml_$i;)
	$(foreach i,$(MP),acol -read mp_$i markov_chain_0 -index -select-rows "N>511" -delete-col N -internal files/mp_$i;)
	$(foreach i,$(NL),acol -read nl_$i markov_chain_0 -index -select-rows "N>511" -delete-col N -internal files/nl_$i;)
	$(foreach i,$(NP),acol -read np_$i markov_chain_0 -index -select-rows "N>511" -delete-col N -internal files/np_$i;)

file-all-con:
	acol -read ml_31 markov_chain_0 $(foreach i,$(ML),-cat ml_$i) -internal files/main_ML.o2
	acol -read mp_32 markov_chain_0 $(foreach i,$(MP),-cat mp_$i) -internal files/main_MP.o2
	acol -read nl_28 markov_chain_0 $(foreach i,$(NL),-cat nl_$i) -internal files/main_NL.o2
	acol -read np_30 markov_chain_0 $(foreach i,$(NP),-cat np_$i) -internal files/main_NP.o2

file-cp:
	$(foreach m,$(models),\
        $(foreach r,$(rank3),\
            $(if $(filter $(m),ml),\
                cp $(m)_59b $(m)_59b$(r);,\
                cp $(m)_50a10 $(m)_50a$(r);\
            )\
        )\
    )

R = 0
T = 5
W = 30

#file-init:
#	acol -read mp_51a$R markov_chain_0 \
	-select-rows "thread==$T && walker==$W" \
	-index -function "(-N)" N2 -sort N2 -index N3 -select-rows "N3<1" \
	-function 0 thread -function 0 walker \
	-delete-col N -delete-col N2 -delete-col N3 -internal mpa$R_$T_$W

file-view:
	acol -read np_50a$R markov_chain_0 \
	-index -select-rows "thread==$T && walker==$W" \
	-select N thread walker log_wgt -o

N = 5284
#N = 1910

file-point:
	acol -read np_50a$R markov_chain_0 \
	-index -select-rows "N==$N" \
	-function 0 thread -function 0 walker \
	-delete-col N -internal npa$R_$T_$W

file-init:
	acol -read np_50a$R markov_chain_0 \
	-index -select-rows "N!=$N" \
	-delete-col N -internal np_50a$R

file-rm:
	rm -f *_out *_scr

clean:
	rm -f figs/*.png figs/*.pdf files/*.o2 
