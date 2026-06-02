
MODEL =
title =

INDICES = $(shell seq 0 99)
HBARC = 197.326980
MN = 939

# Base file - create either by concatenating, or copying, or processing raw files
file-main:
	acol -read $(MODEL)_all markov_chain_0 -select-rows "mult>0.5" \
	-internal tmp/tmp_$(MODEL).o2
	acol -h5-copy tmp/tmp_$(MODEL).o2 tmp/main_$(MODEL).o2
	rm -f tmp/tmp_$(MODEL).o2

# Requires main_$(MODEL).o2
file-main2:
	acol -read tmp/main_$(MODEL).o2 markov_chain_0 \
	-function "10^log10_width_NS" width_NS \
	-function "10^log10_width_WD" width_WD \
	-function "10^log10_width_LMS" width_LM \
	-function "mean_LMS" mean_LM -delete-col mean_LMS \
	-function "skewness_LMS" skewness_LM \
	-delete-col skewness_LMS \
	-function "1/sqrt(2*3.1416)/width_NS" cf \
	-function "exp(-0.5*((1.37-mean_NS)/width_NS)^2)" pdf \
	-function "1+erf((1.37-mean_NS)*(skewness_NS)/width_NS/sqrt(2))" cdf \
	-function "cf*pdf*cdf" SN_NS \
	-function "1/sqrt(2*3.1416)/width_WD" cf \
	-function "exp(-0.5*((1.55-mean_WD)/width_WD)^2)" pdf \
	-function "1+erf((1.55-mean_WD)*(skewness_WD)/width_WD/sqrt(2))" cdf \
	-function "cf*pdf*cdf" SN_WD \
	-function "1/sqrt(2*3.1416)/width_LM" cf \
	-function "exp(-0.5*((1.62-mean_LM)/width_LM)^2)" pdf \
	-function "1+erf((1.62-mean_LM)*(skewness_LM)/width_LM/sqrt(2))" cdf \
	-function "cf*pdf*cdf" SN_LM \
	-function "log(dpdM)" log_dpdM \
	-function "P_23*$(HBARC)" Phc_23 \
	-select mult ":mean_*" ":width_*" ":skewness_*" ":SN_*" \
	Phc_23 R_42 R_43 R_48 R_51 R_max M_max \
	-internal tmp/main2_$(MODEL).o2

XLO = 0.8
XHI = 3.0
YLO = 0.08
YHI = 4.0
DX = 99.0

# Requires main2_$(MODEL).o2
file-dist:
	acol -read tmp/main2_$(MODEL).o2 markov_chain_0 \
	$(foreach i,$(INDICES),-function "1/sqrt(2*3.1416)/width_NS" cf \
	-function "exp(-0.5*(($i*0.03-mean_NS)/width_NS)^2)" pdf \
	-function "1+erf(($i*0.03-mean_NS)*(skewness_NS)/width_NS/sqrt(2))" cdf \
	-function "if($i*0.03<=M_max,cf*pdf*cdf,0.0)" SN_NS_$i) \
	$(foreach i,$(INDICES),-function "1/sqrt(2*3.1416)/width_WD" cf \
	-function "exp(-0.5*(($i*0.03-mean_WD)/width_WD)^2)" pdf \
	-function "1+erf(($i*0.03-mean_WD)*(skewness_WD)/width_WD/sqrt(2))" cdf \
	-function "if($i*0.03<=M_max,cf*pdf*cdf,0.0)" SN_WD_$i) \
	$(foreach i,$(INDICES),-function "1/sqrt(2*3.1416)/width_LM" cf \
	-function "exp(-0.5*(($i*0.03-mean_LM)/width_LM)^2)" pdf \
	-function "1+erf(($i*0.03-mean_LM)*(skewness_LM)/width_LM/sqrt(2))" cdf \
	-function "if($i*0.03<=M_max,cf*pdf*cdf,0.0)" SN_LM_$i) \
	-delete-col SN_NS -delete-col SN_WD -delete-col SN_LM \
	-select ":SN_*" -set obj_name gm_sn -internal tmp/dist_$(MODEL).o2
	acol -read tmp/dist_$(MODEL).o2 gm_sn -ser-hist-t3d \
	"grid:$(XLO),$(XHI),($(XHI)-$(XLO))/$(DX)" x gm \
	"grid:$(YLO),$(YHI),($(YHI)-$(YLO))/$(DX)" auto sn \
	"SN_NS_*" prob_ns -type -list \
	-set obj_name sn_ns -internal tmp/dist_$(MODEL).o2 \
	-to-hist-2d prob_ns \
	-contours frac 0.68 tmp/dist_$(MODEL).o2 c68_ns \
	-contours frac 0.95 tmp/dist_$(MODEL).o2 c95_ns \
	-read tmp/dist_$(MODEL).o2 gm_sn -ser-hist-t3d \
	"grid:$(XLO),$(XHI),($(XHI)-$(XLO))/$(DX)" x gm \
	"grid:$(YLO),$(YHI),($(YHI)-$(YLO))/$(DX)" auto sn \
	"SN_WD_*" prob_wd -type -list \
	-set obj_name sn_wd -internal tmp/dist_$(MODEL).o2 \
	-to-hist-2d prob_wd \
	-contours frac 0.68 tmp/dist_$(MODEL).o2 c68_wd \
	-contours frac 0.95 tmp/dist_$(MODEL).o2 c95_wd \
	-read tmp/dist_$(MODEL).o2 gm_sn -ser-hist-t3d \
	"grid:$(XLO),$(XHI),($(XHI)-$(XLO))/$(DX)" x gm \
	"grid:$(YLO),$(YHI),($(YHI)-$(YLO))/$(DX)" auto sn \
	"SN_LM_*" prob_lm -type -list \
	-set obj_name sn_lm -internal tmp/dist_$(MODEL).o2 \
	-to-hist-2d prob_lm \
	-contours frac 0.68 tmp/dist_$(MODEL).o2 c68_lm \
	-contours frac 0.95 tmp/dist_$(MODEL).o2 c95_lm  

file-dist2:
	acol -read tmp/main2_$(MODEL).o2 markov_chain_0 \
	$(foreach i,$(INDICES),-function "1/sqrt(2*3.1416)/width_NS" cf \
	-function "exp(-0.5*(($i*0.03-mean_NS)/width_NS)^2)" pdf \
	-function "1+erf(($i*0.03-mean_NS)*(skewness_NS)/width_NS/sqrt(2))" cdf \
	-function "if(($i*0.03<=M_max) && ($i*0.03>2.0),cf*pdf*cdf,0.0)" SN_NS_$i) \
	$(foreach i,$(INDICES),-function "1/sqrt(2*3.1416)/width_WD" cf \
	-function "exp(-0.5*(($i*0.03-mean_WD)/width_WD)^2)" pdf \
	-function "1+erf(($i*0.03-mean_WD)*(skewness_WD)/width_WD/sqrt(2))" cdf \
	-function "if(($i*0.03<=M_max) && ($i*0.03>2.0),cf*pdf*cdf,0.0)" SN_WD_$i) \
	$(foreach i,$(INDICES),-function "1/sqrt(2*3.1416)/width_LM" cf \
	-function "exp(-0.5*(($i*0.03-mean_LM)/width_LM)^2)" pdf \
	-function "1+erf(($i*0.03-mean_LM)*(skewness_LM)/width_LM/sqrt(2))" cdf \
	-function "if(($i*0.03<=M_max) && ($i*0.03>2.0),cf*pdf*cdf,0.0)" SN_LM_$i) \
	-delete-col SN_NS -delete-col SN_WD -delete-col SN_LM \
	-select ":SN_*" M_max -set obj_name gm_sn -internal tmp/dist2_$(MODEL).o2

file-esmr:
	acol -read tmp/main_$(MODEL).o2 markov_chain_0 \
	$(foreach i,$(INDICES),-function P_$i*$(HBARC) Phc_$i) \
	-set obj_name ed_pr -select ":Phc_*" -internal tmp/esmr_$(MODEL).o2 \
	-read tmp/esmr_$(MODEL).o2 ed_pr -ser-hist-t3d \
	"grid:0.3*$(HBARC),10.0*$(HBARC),((10.0-0.3)*$(HBARC))/99.0" x ed \
	"grid:1.0,2.0e3,(2.0e3)^(1.0/99.0),log" auto pr \
	"Phc_*" prob -type -list \
	-set obj_name ed_pr_3d -internal tmp/esmr_$(MODEL).o2 \
	-to-hist-2d prob \
	-contours frac 0.68 tmp/esmr_$(MODEL).o2 c68_es \
	-contours frac 0.95 tmp/esmr_$(MODEL).o2 c95_es \
	-read tmp/main_$(MODEL).o2 markov_chain_0 -delete-col R_max \
	-set obj_name gm_r -select ":R_*" -internal tmp/esmr_$(MODEL).o2 \
	-read tmp/esmr_$(MODEL).o2 gm_r -ser-hist-t3d \
	"grid:0.2,3.0,(3.0-0.2)/99.0" y gm \
	"grid:9.0,15.0,(15.0-9.0)/99.0" auto r \
	"R_*" prob -type -list \
	-set obj_name gm_r_3d -internal tmp/esmr_$(MODEL).o2 \
	-to-hist-2d prob \
	-contours frac 0.68 tmp/esmr_$(MODEL).o2 c68_mr \
	-contours frac 0.95 tmp/esmr_$(MODEL).o2 c95_mr

file-eoa:
	acol -read tmp/main_$(MODEL).o2 markov_chain_0 \
	$(foreach i,$(INDICES),-function "EoA_$i*$(HBARC)" EoA_hc_$i) \
	-select ":EoA_hc_*" -set obj_name nb_eoa -internal tmp/eoa_$(MODEL).o2 \
	-read tmp/eoa_$(MODEL).o2 nb_eoa -ser-hist-t3d \
	"grid:0.04,1.24,(1.24-0.04)/99.0" x nb \
	"grid:1.0,4.32*$(HBARC),(4.32*$(HBARC))^(1/99),log" auto ea \
	"EoA_hc_*" prob -type -list \
	-set obj_name nb_eoa_3d -internal tmp/eoa_$(MODEL).o2 \
	-to-hist-2d prob \
	-contours frac 0.68 tmp/eoa_$(MODEL).o2 c68 \
	-contours frac 0.95 tmp/eoa_$(MODEL).o2 c95

file-pnb:
	acol -read tmp/main_$(MODEL).o2 markov_chain_0 \
	$(foreach i,$(INDICES),-function Pnb_$i*$(HBARC) Pnb_hc_$i) \
	-select ":Pnb_hc_*" -set obj_name nb_pr -internal tmp/pnb_$(MODEL).o2 \
	-read tmp/pnb_$(MODEL).o2 nb_pr -ser-hist-t3d \
	"grid:0.04,1.24,(1.24-0.04)/99.0" x nb \
	"grid:1.0,2.0e3,(2.0e3)^(1/99),log" auto pr \
	"Pnb_hc_*" prob -type -list \
	-set obj_name nb_pr_3d -internal tmp/pnb_$(MODEL).o2 \
	-to-hist-2d prob \
	-contours frac 0.68 tmp/pnb_$(MODEL).o2 c68 \
	-contours frac 0.95 tmp/pnb_$(MODEL).o2 c95

file-cs2:
	acol -read tmp/main_$(MODEL).o2 markov_chain_0 \
	$(foreach i,$(INDICES),-function "if(cs2_$i>1.0,1.0,cs2_$i)" cs2_n_$i) \
	-set obj_name ed_cs2 -select ":cs2_n_*" -internal tmp/cs2_$(MODEL).o2 \
	-read tmp/cs2_$(MODEL).o2 ed_cs2 -ser-hist-t3d \
	"grid:0.3*$(HBARC),1550.0,(1550.0-0.3*$(HBARC))/99.0" x ed \
	"grid:0.0,1.0,(1.0-0.0)/99.0" auto cs2 \
	"cs2_*" prob -type -list \
	-set obj_name ed_cs2_3d -internal tmp/cs2_$(MODEL).o2 \
	-to-hist-2d prob \
	-contours frac 0.68 tmp/cs2_$(MODEL).o2 c68 \
	-contours frac 0.95 tmp/cs2_$(MODEL).o2 c95

file-lambdat:
	acol -read tmp/main_$(MODEL).o2 markov_chain_0 \
	-select ":Lambda_bar_*" -set obj_name m_lt -internal tmp/lambdat_$(MODEL).o2 \
	-read tmp/lambdat_$(MODEL).o2 m_lt -ser-hist-t3d \
	"grid:0.2,3.0,(3.0-0.2)/99.0" x gm \
	"grid:1,5000,(5000-1)^(1/99),log" auto lt \
	"Lambda_bar_*" prob -type -list \
	-set obj_name m_lt_3d -internal tmp/lambdat_$(MODEL).o2 \
	-to-hist-2d prob \
	-contours frac 0.68 tmp/lambdat_$(MODEL).o2 c68 \
	-contours frac 0.90 tmp/lambdat_$(MODEL).o2 c90 \
	-contours frac 0.95 tmp/lambdat_$(MODEL).o2 c95

file-lambdat-all:
	make MODEL="ml" file-lambdat && make MODEL="mp" file-lambdat && make MODEL="nl" file-lambdat && make MODEL="np" file-lambdat

file-eoa-all:
	make MODEL="ml" file-eoa && make MODEL="mp" file-eoa && make MODEL="nl" file-eoa && make MODEL="np" file-eoa

file-pnb-all:
	make MODEL="ml" file-pnb && make MODEL="mp" file-pnb && make MODEL="nl" file-pnb && make MODEL="np" file-pnb

file-cs2-all:
	make MODEL=ml file-cs2 && make MODEL=mp file-cs2 && make MODEL=nl file-cs2 && make MODEL=np file-cs2
	
file-esmr-all:
	make MODEL="ml" file-esmr && make MODEL="mp" file-esmr && make MODEL="nl" file-esmr && make MODEL="np" file-esmr

file-main2-all:
	make MODEL="ml" file-main2 && make MODEL="mp" file-main2 && make MODEL="nl" file-main2 && make MODEL="np" file-main2

file-dist-all:
	make MODEL="ml" file-dist && make MODEL="mp" file-dist && make MODEL="nl" file-dist && make MODEL="np" file-dist

file-dist2-all:
	make MODEL="ml" file-dist2 && make MODEL="mp" file-dist2 && make MODEL="nl" file-dist2 && make MODEL="np" file-dist2

file-all:
	make file-main2-all && make file-dist-all && make file-esmr-all && make file-eoa-all && make file-pnb-all && make file-cs2-all && make file-lambdat-all

RANKS := $(shell seq 0 15)
MODEL := np
INDEX := 88

.PHONY: remove-zeros
remove-zeros:
	@$(foreach r,$(RANKS), \
	  echo "Removing zero mults in $(MODEL)_$(r)_out"; \
	  cmd='acol -read $(MODEL)_$(r)_out markov_chain_0 \
	       -select-rows "mult>0" \
	       -internal $(MODEL)_$(r)_out'; \
	  echo "$$cmd"; \
	  eval $$cmd || echo "Failed on $(MODEL)_$(r)_out"; \
	)

.PHONY: copy-hdf5
copy-hdf5:
	@$(foreach r,$(RANKS), \
	  echo "Copying hdf5 file $(MODEL)_$(r)_out"; \
	  cmd="acol -h5-copy $(MODEL)_$(r)_out \
	       $(MODEL)_$(INDEX)r$(r)"; \
	  echo "$$cmd"; \
	  $$cmd || echo "Failed on $(MODEL)_$(r)_out"; \
	)

RANKS1 := $(shell seq 1 15)

.PHONY: remove-init
remove-init:
	@$(foreach r,$(RANKS), \
	  echo "Removing initial rows (N<512) in files/$(MODEL)_$(INDEX)r$(r)"; \
	  cmd='acol -read files/$(MODEL)_$(INDEX)r$(r) markov_chain_0 \
	       -index -select-rows "N>=512" -delete-col N \
	       -internal thin/$(MODEL)_$(INDEX)r$(r)'; \
	  echo "$$cmd"; \
	  eval $$cmd || echo "Failed on thin/$(MODEL)_$(INDEX)r$(r)"; \
	)

.PHONY: cat-ranks
cat-ranks:
	@echo "Concatenating all ranks into $(MODEL)_$(INDEX)"
	@cmd="acol -read thin/$(MODEL)_$(INDEX)r0 markov_chain_0 $(foreach j,$(RANKS1),-cat thin/$(MODEL)_$(INDEX)r$(j)) -internal thin/$(MODEL)_$(INDEX)"; \
	echo "$$cmd"; \
	eval $$cmd || echo "Cat failed for $(MODEL)_$(INDEX)"

MODELS := ml mp nl np

.PHONY: thin-to-tmp
thin-to-tmp:
	@echo "Copying thinned files to tmp/"
	cp -v thin/ml_thin tmp/main_ml.o2
	cp -v thin/mp_thin tmp/main_mp.o2
	cp -v thin/nl_thin tmp/main_nl.o2
	cp -v thin/np_thin tmp/main_np.o2

.PHONY: file-pars
file-pars:
	@echo "Extracting parameters from thinned models"
	@for m in $(MODELS); do \
	  if [ "$$m" = ml ] || [ "$$m" = nl ]; then \
	    EOSP="csq1 trans1 csq2 trans2 csq3"; \
	  else \
	    EOSP="exp1 trans1 exp2 trans2 exp3"; \
	  fi; \
	  IN="thin/$${m}_thin"; \
	  OUT="thin/$${m}_thin2"; \
	  cmd="acol -read $$IN markov_chain_0 -select a alpha param_S param_L $$EOSP M_chirp_det q z_cdf m1_gw19 mf_6304 mf_6397 mf_M13 mf_M28 mf_M30 mf_wCen mf_X7 mf_1810b mf_1724b mf_1702 mf_0030 mf_0740 mean_NS log10_width_NS skewness_NS mean_WD log10_width_WD skewness_WD mean_LMS log10_width_LMS skewness_LMS M_J0453p M_J0453c M_J1906p M_J1906c M_B1534p M_B1534c M_B1913p M_B1913c M_B2127p M_B2127c M_J0737A M_J0737B M_J1756p M_J1756c M_J1807p M_J1807c M_J1518p M_J1518c M_J1811p M_J1811c M_J1829p M_J1829c M_J2045 M_J2053 M_J1713 M_B1855 M_J0751 M_J1141 M_J1738 M_J1614 M_J0348 M_J2222 M_J2234 M_J1949 M_J1012 M_J0437 M_J1909 M_J1802 M_J1911 M_J2043 M_J0337 M_J1946 M_J1918 M_J1600 M_J0621 M_B2303 M_J0024 M_J0514 M_B1516 M_J1748I M_J1748J M_B1802 M_B1911 M_J0740 M_CygX2 M_XTEJ2123 M_4U1822 M_HerX1 M_2S0921 -internal $$OUT"; \
	  echo "$$cmd"; \
	  eval $$cmd || echo "Failed on $$m"; \
	done

.PHONY: clean-thin
clean-thin:
	@echo "Removing files: thin/$(MODEL)_$(INDEX)r* and thin/$(MODEL)_thin*"
	rm -fv thin/$(MODEL)_$(INDEX)r* thin/$(MODEL)_thin thin/$(MODEL)_thin2

.PHONY: clean-tmp
clean-tmp:
	@echo "Removing files: tmp/*.o2"
	rm -fv tmp/*.o2

.PHONY:
all:
	make remove-zeros && make copy-hdf5 && cp -v $(MODEL)_$(INDEX)r* files/ && \
	make remove-init && make cat-ranks && rm -fv thin/$(MODEL)_$(INDEX)r*