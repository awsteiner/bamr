
model =

title =

INDEX = 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 \
	20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 \
	40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 \
	60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 \
	80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99

INDEX2 = 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 \
	76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100 101 102 \
	103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 \
	124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 144 \
	145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 162 163 164 165 \
	166 167 168 169 170 171 172 173 174 175 176 177 178 179 180 181 182 183 184 185 186 \
	187 188 189 190 191 192 193 194 195 196 197 198 199 200 201 202 203 204 205 206 207 \
	208 209 210 211 212 213 214 215 216 217 218 219 220 221 222 223 224 225 226 227 228 \
	229 230 231 232 233 234 235 236 237 238 239 240 241 242 243 244 245 246 247 248 249 \

THREAD = 0 

WALKER = 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 \
	20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 \
	40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 \
	60 61 62 63

IWALK = 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 \
		31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 \
		61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 \
		91 92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 \
		116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 \
		139 140 141 142 143 144 145 146 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161 \
		162 163 164 165 166 167 168 169 170 171 172 173 174 175 176 177 178 179 180 181 182 183 184 \
		185 186 187 188 189 190 191 192 193 194 195 196 197 198 199 200 201 202 203 204 205 206 207 \
		208 209 210 211 212 213 214 215 216 217 218 219 220 221 222 223 224 225 226 227 228 229 230 \
		231 232 233 234 235 236 237 238 239 240 241 242 243 244 245 246 247 248 249 250 251 252 253 \
		254 255 256 257 258 259 260 261 262 263 264 265 266 267 268 269 270 271 272 273 274 275 276 \
		277 278 279 280 281 282 283 284 285 286 287 288 289 290 291 292 293 294 295 296 297 298 299 \
		300 301 302 303 304 305 306 307 308 309 310 311 312 313 314 315 316 317 318 319 320 321 322 \
		323 324 325 326 327 328 329 330 331 332 333 334 335 336 337 338 339 340 341 342 343 344 345 \
		346 347 348 349 350 351 352 353 354 355 356 357 358 359 360 361 362 363 364 365 366 367 368 \
		369 370 371 372 373 374 375 376 377 378 379 380 381 382 383 384 385 386 387 388 389 390 391 \
		392 393 394 395 396 397 398 399 400 401 402 403 404 405 406 407 408 409 410 411 412 413 414 \
		415 416 417 418 419 420 421 422 423 424 425 426 427 428 429 430 431 432 433 434 435 436 437 \
		438 439 440 441 442 443 444 445 446 447 448 449 450 451 452 453 454 455 456 457 458 459 460 \
		461 462 463 464 465 466 467 468 469 470 471 472 473 474 475 476 477 478 479 480 481 482 483 \
		484 485 486 487 488 489 490 491 492 493 494 495 496 497 498 499 500 501 502 503 504 505 506 \
		507 508 509 510

new-file:
	acol -read mp2_31 markov_chain_0 $(foreach i, $(IWALK),-cat mp2_31) -internal mp3_31

# Base file - create either by concatenating, or copying, or processing raw files
file-main:
	acol -read out_$(model) markov_chain_0 -select-rows "mult>0.5" \
	-internal tmp/tmp_$(model).o2
	acol -h5-copy tmp/tmp_$(model).o2 tmp/main_$(model).o2
	rm -f tmp/tmp_$(model).o2

# Requires main_$(model).o2
file-esmr:
	acol -read tmp/main_$(model).o2 markov_chain_0 \
	$(foreach i,$(INDEX),-function P_$i*197.32698 Phc_$i) \
	-set obj_name ed_pr -select ":Phc_*" -internal tmp/esmr_$(model).o2 \
	-read tmp/esmr_$(model).o2 ed_pr -ser-hist-t3d \
	"grid:0.3*197.32698,10.0*197.32698,((10.0-0.3)*197.32698)/99.0" x ed \
	"grid:1.0,2.0e3,(2.0e3)^(1.0/99.0),log" auto pr \
	"Phc_*" prob -type -list \
	-set obj_name ed_pr_3d -internal tmp/esmr_$(model).o2 \
	-to-hist-2d prob \
	-contours frac 0.68 tmp/esmr_$(model).o2 c68_es \
	-contours frac 0.95 tmp/esmr_$(model).o2 c95_es \
	-read tmp/main_$(model).o2 markov_chain_0 -delete-col R_max \
	-set obj_name gm_r -select ":R_*" -internal tmp/esmr_$(model).o2 \
	-read tmp/esmr_$(model).o2 gm_r -ser-hist-t3d \
	"grid:0.2,3.0,(3.0-0.2)/99.0" y gm \
	"grid:9.0,15.0,(15.0-9.0)/99.0" auto r \
	"R_*" prob -type -list \
	-set obj_name gm_r_3d -internal tmp/esmr_$(model).o2 \
	-to-hist-2d prob \
	-contours frac 0.68 tmp/esmr_$(model).o2 c68_mr \
	-contours frac 0.95 tmp/esmr_$(model).o2 c95_mr

# Requires main_$(model).o2
file-main2:
	acol -read tmp/main_$(model).o2 markov_chain_0 \
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
	-internal tmp/main2_$(model).o2

XLO = 0.0
XHI = 3.0
YLO = 0.04
YHI = 4.0
DX = 99.0

# Requires main2_$(model).o2
file-dist:
	acol -read tmp/main2_$(model).o2 markov_chain_0 \
	$(foreach i,$(INDEX),-function "1/sqrt(2*3.1416)/width_NS" cf \
	-function "exp(-0.5*(($i*0.03-mean_NS)/width_NS)^2)" pdf \
	-function "1+erf(($i*0.03-mean_NS)*(skewness_NS)/width_NS/sqrt(2))" cdf \
	-function "cf*pdf*cdf" SN_NS_$i) \
	$(foreach i,$(INDEX),-function "1/sqrt(2*3.1416)/width_WD" cf \
	-function "exp(-0.5*(($i*0.03-mean_WD)/width_WD)^2)" pdf \
	-function "1+erf(($i*0.03-mean_WD)*(skewness_WD)/width_WD/sqrt(2))" cdf \
	-function "cf*pdf*cdf" SN_WD_$i) \
	$(foreach i,$(INDEX),-function "1/sqrt(2*3.1416)/width_LM" cf \
	-function "exp(-0.5*(($i*0.03-mean_LM)/width_LM)^2)" pdf \
	-function "1+erf(($i*0.03-mean_LM)*(skewness_LM)/width_LM/sqrt(2))" cdf \
	-function "cf*pdf*cdf" SN_LM_$i) \
	-select ":SN_*" -set obj_name gm_sn -internal tmp/dist_$(model).o2
	acol -read tmp/dist_$(model).o2 gm_sn -ser-hist-t3d \
	"grid:$(XLO),$(XHI),($(XHI)-$(XLO))/$(DX)" x gm \
	"grid:$(YLO),$(YHI),($(YHI)-$(YLO))/$(DX)" auto sn \
	"SN_NS_*" prob_ns -type -list \
	-set obj_name sn_ns -internal tmp/dist_$(model).o2 \
	-to-hist-2d prob_ns \
	-contours frac 0.68 tmp/dist_$(model).o2 c68_ns \
	-contours frac 0.95 tmp/dist_$(model).o2 c95_ns \
	-read tmp/dist_$(model).o2 gm_sn -ser-hist-t3d \
	"grid:$(XLO),$(XHI),($(XHI)-$(XLO))/$(DX)" x gm \
	"grid:$(YLO),$(YHI),($(YHI)-$(YLO))/$(DX)" auto sn \
	"SN_WD_*" prob_wd -type -list \
	-set obj_name sn_wd -internal tmp/dist_$(model).o2 \
	-to-hist-2d prob_wd \
	-contours frac 0.68 tmp/dist_$(model).o2 c68_wd \
	-contours frac 0.95 tmp/dist_$(model).o2 c95_wd \
	-read tmp/dist_$(model).o2 gm_sn -ser-hist-t3d \
	"grid:$(XLO),$(XHI),($(XHI)-$(XLO))/$(DX)" x gm \
	"grid:$(YLO),$(YHI),($(YHI)-$(YLO))/$(DX)" auto sn \
	"SN_LM_*" prob_lm -type -list \
	-set obj_name sn_lm -internal tmp/dist_$(model).o2 \
	-to-hist-2d prob_lm \
	-contours frac 0.68 tmp/dist_$(model).o2 c68_lm \
	-contours frac 0.95 tmp/dist_$(model).o2 c95_lm 

splot-eos1:
	o2graph -subplots 2 1 sharex=True \
	-subadj hspace=0.0,left=0.14,right=0.84,top=0.98,bottom=0.12 \
	-selax 0 -set logy 1 -read tmp/esmr_NL.o2 ed_pr_3d \
	-den-plot prob pcm=True,cmap=pink_r \
	-read tmp/esmr_NL.o2 c68_es -plot lw=1,color=purple \
	-read tmp/esmr_NL.o2 c95_es -plot lw=1,color=orange \
	-ttext 0.1 0.92 "NL" \
	-selax 1 -set logy 1 -read tmp/esmr_NP.o2 ed_pr_3d \
	-den-plot prob pcm=True,cmap=pink_r \
	-read tmp/esmr_NP.o2 c68_es -plot lw=1,color=purple \
	-read tmp/esmr_NP.o2 c95_es -plot lw=1,color=orange \
	-xtitle "Energy Density, $$ \epsilon~[\mathrm{MeV}/\mathrm{fm}^3] $$" \
	-ytitle ".$$ \qquad \qquad \qquad \qquad \qquad $$ Pressure, $$ P~[\mathrm{MeV}/\mathrm{fm}^3] $$" \
	-ttext 0.1 0.92 "NP" \
	-addcbar 0.86 0.12 0.03 0.86 -selax \
	-save figs/eos_s1.png -show &

splot-eos2:
	o2graph -subplots 2 1 sharex=True \
	-subadj hspace=0.0,left=0.14,right=0.84,top=0.98,bottom=0.12 \
	-selax 0 -set logy 1 -read tmp/esmr_ML.o2 ed_pr_3d \
	-den-plot prob pcm=True,cmap=pink_r \
	-read tmp/esmr_ML.o2 c68_es -plot lw=1,color=purple \
	-read tmp/esmr_ML.o2 c95_es -plot lw=1,color=orange \
	-ttext 0.1 0.92 "ML" \
	-selax 1 -set logy 1 -read tmp/esmr_MP.o2 ed_pr_3d \
	-den-plot prob pcm=True,cmap=pink_r \
	-read tmp/esmr_MP.o2 c68_es -plot lw=1,color=purple \
	-read tmp/esmr_MP.o2 c95_es -plot lw=1,color=orange \
	-xtitle "Energy Density, $$ \epsilon~[\mathrm{MeV}/\mathrm{fm}^3] $$" \
	-ytitle ".$$ \qquad \qquad \qquad \qquad \qquad $$ Pressure, $$ P~[\mathrm{MeV}/\mathrm{fm}^3] $$" \
	-ttext 0.1 0.92 "MP" \
	-addcbar 0.86 0.12 0.03 0.86 -selax \
	-save figs/eos_s2.png -show &

splot-eos:
	o2graph -subplots 4 1 sharex=True \
	-subadj hspace=0.0,left=0.16,right=0.8,top=0.99,bottom=0.05 \
	-selax 0 -set logy 1 -read tmp/esmr_NL.o2 ed_pr_3d \
	-den-plot prob pcm=True,cmap=pink_r \
	-read tmp/esmr_NL.o2 c68_es -plot lw=1,color=purple \
	-read tmp/esmr_NL.o2 c95_es -plot lw=1,color=orange \
	-ttext 0.1 0.92 "NL" \
	-selax 1 -set logy 1 -read tmp/esmr_NP.o2 ed_pr_3d \
	-den-plot prob pcm=True,cmap=pink_r \
	-read tmp/esmr_NP.o2 c68_es -plot lw=1,color=purple \
	-read tmp/esmr_NP.o2 c95_es -plot lw=1,color=orange \
	-ttext 0.1 0.92 "NP" \
	-selax 2 -set logy 1 -read tmp/esmr_ML.o2 ed_pr_3d \
	-den-plot prob pcm=True,cmap=pink_r \
	-read tmp/esmr_ML.o2 c68_es -plot lw=1,color=purple \
	-read tmp/esmr_ML.o2 c95_es -plot lw=1,color=orange \
	-ttext 0.1 0.92 "ML" \
	-ytitle ".$$ \qquad \qquad \qquad \qquad \qquad \qquad \quad $$ Pressure, $$ P~[\mathrm{MeV}/\mathrm{fm}^3] $$" \
	-selax 3 -set logy 1 -read tmp/esmr_MP.o2 ed_pr_3d \
	-den-plot prob pcm=True,cmap=pink_r \
	-read tmp/esmr_MP.o2 c68_es -plot lw=1,color=purple \
	-read tmp/esmr_MP.o2 c95_es -plot lw=1,color=orange \
	-ttext 0.1 0.92 "MP" \
	-xtitle "Energy Density, $$ \epsilon~[\mathrm{MeV}/\mathrm{fm}^3] $$" \
	-addcbar 0.82 0.05 0.03 0.94 -selax \
	-show &

splot-mvsr1:
	o2graph -subplots 2 1 sharex=True \
	-subadj hspace=0.0,left=0.14,right=0.84,top=0.98,bottom=0.12 \
	-selax 0 -read tmp/esmr_NL.o2 gm_r_3d -den-plot prob cmap=pink_r \
	-read tmp/esmr_NL.o2 c68_mr -plot lw=1,color=purple \
	-read tmp/esmr_NL.o2 c95_mr -plot lw=1,color=orange \
	-ttext 0.1 0.92 "NL" \
	-selax 1 -read tmp/esmr_NP.o2 gm_r_3d -den-plot prob cmap=pink_r \
	-read tmp/esmr_NP.o2 c68_mr -plot lw=1,color=purple \
	-read tmp/esmr_NP.o2 c95_mr -plot lw=1,color=orange \
	-ttext 0.1 0.92 "NP" \
	-xtitle "Radius, $$ R~[\mathrm{km}] $$" \
	-ytitle ".$$ \qquad \qquad \qquad \qquad \qquad \quad $$ Mass, $$ M~[\mathrm{M}_\odot] $$" \
	-addcbar 0.86 0.12 0.03 0.86 -selax \
	-save figs/mvsr_s1.png -show &

splot-mvsr2:
	o2graph -subplots 2 1 sharex=True \
	-subadj hspace=0.0,left=0.14,right=0.84,top=0.98,bottom=0.12 \
	-selax 0 -read tmp/esmr_ML.o2 gm_r_3d -den-plot prob cmap=pink_r \
	-read tmp/esmr_ML.o2 c68_mr -plot lw=1,color=purple \
	-read tmp/esmr_ML.o2 c95_mr -plot lw=1,color=orange \
	-ttext 0.1 0.92 "ML" \
	-selax 1 -read tmp/esmr_MP.o2 gm_r_3d -den-plot prob cmap=pink_r \
	-read tmp/esmr_MP.o2 c68_mr -plot lw=1,color=purple \
	-read tmp/esmr_MP.o2 c95_mr -plot lw=1,color=orange \
	-xtitle "Radius, $$ R~[\mathrm{km}] $$" \
	-ytitle ".$$ \qquad \qquad \qquad \qquad \qquad \quad $$ Mass, $$ M~[\mathrm{M}_\odot] $$" \
	-ttext 0.1 0.92 "MP" \
	-addcbar 0.86 0.12 0.03 0.86 -selax \
	-save figs/mvsr_s2.png -show &

splot-mvsr:
	o2graph -subplots 4 1 sharex=True \
	-subadj hspace=0.0,left=0.16,right=0.8,top=0.99,bottom=0.05 \
	-selax 0 -read tmp/esmr_NL.o2 gm_r_3d \
	-den-plot prob cmap=pink_r \
	-read tmp/esmr_NL.o2 c68_mr -plot lw=1,color=purple \
	-read tmp/esmr_NL.o2 c95_mr -plot lw=1,color=orange \
	-ttext 0.1 0.92 "NL" \
	-selax 1 -read tmp/esmr_NP.o2 gm_r_3d \
	-den-plot prob cmap=pink_r \
	-read tmp/esmr_NP.o2 c68_mr -plot lw=1,color=purple \
	-read tmp/esmr_NP.o2 c95_mr -plot lw=1,color=orange \
	-ttext 0.1 0.92 "NP" \
	-selax 2 -read tmp/esmr_ML.o2 gm_r_3d \
	-den-plot prob cmap=pink_r \
	-read tmp/esmr_ML.o2 c68_mr -plot lw=1,color=purple \
	-read tmp/esmr_ML.o2 c95_mr -plot lw=1,color=orange \
	-ttext 0.1 0.92 "ML" \
	-ytitle ".$$ \qquad \qquad \qquad \qquad \qquad \qquad \quad $$ Mass, $$ M~[\mathrm{M}_\odot] $$" \
	-selax 3 -read tmp/esmr_MP.o2 gm_r_3d \
	-den-plot prob cmap=pink_r \
	-read tmp/esmr_MP.o2 c68_mr -plot lw=1,color=purple \
	-read tmp/esmr_MP.o2 c95_mr -plot lw=1,color=orange \
	-ttext 0.1 0.92 "MP" \
	-xtitle "Radius, $$ R~[\mathrm{km}] $$" \
	-addcbar 0.82 0.05 0.03 0.94 -selax \
	-show &

plot-hist-mmax:
	o2graph -read tmp/main2_$(model).o2 markov_chain_0 \
	-to-hist M_max 40 log_wgt -plot density=True \
	-xtitle "M_max [$$ M_\odot $$]" -ytitle "Frequency" \
	-ttext 0.5 0.96 "$(model): $(title)" \
	-save figs/hist_mmax_$(model).png -show &

plot-scat-r-sn:
	o2graph -read tmp/main2_$(model).o2 markov_chain_0 -function 1 w -set colbar 1 \
	-scatter R_43 SN_NS w log_wgt cmap=cividis \
	-xtitle "Radius $$ R_{1.4M_\odot} $$ [km]" \
	-ytitle "Mass Distribution, $$ \mathcal{P}(M=1.4 M_\odot $$)" \
	-ttext 0.5 0.96 "$(model): DNS" -save figs/scat_sn-r_ns_$(model).png -show &
	o2graph -read tmp/main2_$(model).o2 markov_chain_0 -function 1 w -set colbar 1 \
	-scatter R_50 SN_WD w log_wgt cmap=cividis \
	-xtitle "Radius $$ R_{1.6M_\odot} $$ [km]" \
	-ytitle "Mass Distribution, $$ \mathcal{P}(M=1.6 M_\odot $$)" \
	-ttext 0.5 0.96 "$(model): NS-WD" -save figs/scat_sn-r_wd_$(model).png -show &
	o2graph -read tmp/main2_$(model).o2 markov_chain_0 -function 1 w -set colbar 1 \
	-scatter R_53 SN_LM w log_wgt cmap=cividis \
	-xtitle "Radius $$ R_{1.7M_\odot} $$ [km]" \
	-ytitle "Mass Distribution, $$ \mathcal{P}(M=1.7 M_\odot $$)" \
	-ttext 0.5 0.96 "$(model): LMXB" -save figs/scat_sn-r_lm_$(model).png -show &

splot-r-sn:
	o2graph -read tmp/main2_$(model).o2 markov_chain_0 -subplots 3 1 \
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
	o2graph -read tmp/main2_NL.o2 markov_chain_0 \
	-kde-2d-plot R_53 SN_LM \
	-ttext 0.1 0.92 "NL" -save figs/scat-r53_NL.pdf -show &
	o2graph -read tmp/main2_NL.o2 markov_chain_0 \
	-kde-2d-plot Phc_23 SN_LM \
	-save figs/scat-p23_NL.pdf -show &
	o2graph -read tmp/main2_NL.o2 markov_chain_0 \
	-kde-2d-plot M_max SN_LM \
	-save figs/scat-mmax_NL.pdf-show &
	o2graph -read tmp/main2_NP.o2 markov_chain_0 \
	-kde-2d-plot R_53 SN_LM \
	-ttext 0.1 0.92 "NP" \
	-ytitle "$$ \mathcal{P}(M=1.7 M_\odot $$)" -show &
	o2graph -read tmp/main2_NP.o2 markov_chain_0 \
	-kde-2d-plot Phc_23 SN_LM -show &
	o2graph -read tmp/main2_NP.o2 markov_chain_0 \
	-kde-2d-plot M_max SN_LM -show &	
	o2graph -read tmp/main2_ML.o2 markov_chain_0 \
	-kde-2d-plot R_53 SN_LM \
	-ttext 0.1 0.92 "ML" -show &
	o2graph -read tmp/main2_ML.o2 markov_chain_0 \
	-kde-2d-plot Phc_23 SN_LM -show &
	o2graph -read tmp/main2_ML.o2 markov_chain_0 \
	-kde-2d-plot M_max SN_LM -show &
	o2graph -read tmp/main2_MP.o2 markov_chain_0 \
	-kde-2d-plot R_53 SN_LM \
	-ttext 0.1 0.92 "MP" \
	-xtitle "$$ R_{1.7M_\odot} $$ [km]" -show &
	o2graph -read tmp/main2_MP.o2 markov_chain_0 \
	-kde-2d-plot Phc_23 SN_LM \
	-xtitle "$$ P~(\epsilon \approx $$ 500 MeV/fm $$ ^3) $$" -show &
	o2graph -read tmp/main2_MP.o2 markov_chain_0 \
	-kde-2d-plot M_max SN_LM \
	-xtitle "$$ M_\mathrm{max} [M_\odot] $$" -show &

plot-scat-p23-sn:
	o2graph -read tmp/main2_$(model).o2 markov_chain_0 -function 1 w -set colbar 1 \
	-scatter Phc_23 SN_NS w log_wgt cmap=cividis \
	-xtitle "Pressure $$ P~(\epsilon \approx $$ 500 MeV/fm $$ ^3) $$" \
	-ytitle "Mass Distribution, $$ \mathcal{P}(M=1.4 M_\odot $$)" \
	-ttext 0.5 0.96 "$(model): DNS" \
	-save figs/scat_sn-p23_$(model).png -show &
	o2graph -read tmp/main2_$(model).o2 markov_chain_0 -function 1 w -set colbar 1 \
	-scatter Phc_23 SN_WD w log_wgt cmap=cividis \
	-xtitle "Pressure $$ P~(\epsilon \approx $$ 500 MeV/fm $$ ^3) $$" \
	-ytitle "Mass Distribution, $$ \mathcal{P}(M=1.6 M_\odot $$)" \
	-ttext 0.5 0.96 "$(model): NS-WD" \
	-save figs/scat_sn-p23_$(model).png -show &
	o2graph -read tmp/main2_$(model).o2 markov_chain_0 -function 1 w -set colbar 1 \
	-scatter Phc_23 SN_LM w log_wgt cmap=cividis \
	-xtitle "Pressure $$ P~(\epsilon \approx $$ 500 MeV/fm $$ ^3) $$" \
	-ytitle "Mass Distribution, $$ \mathcal{P}(M=1.7 M_\odot $$)" \
	-ttext 0.5 0.96 "$(model): LMXB" \
	-save figs/scat_sn-p23_$(model).png -show &

splot-p23-sn:
	o2graph -read tmp/main2_$(model).o2 markov_chain_0 -subplots 3 1 sharex=True \
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
	o2graph -read tmp/main2_$(model).o2 markov_chain_0 -function 1 w -set colbar 1 \
	-scatter M_max SN_NS w log_wgt cmap=cividis \
	-xtitle "Maximum Mass $$ M_\mathrm{max} [M_\odot] $$" \
	-ytitle "Mass Distribution, $$ \mathcal{P}(M=1.4 M_\odot $$)" \
	-ttext 0.5 0.96 "$(model): DNS" \
	-save figs/scat_sn-mmax_$(model).png -show &
	o2graph -read tmp/main2_$(model).o2 markov_chain_0 -function 1 w -set colbar 1 \
	-scatter M_max SN_WD w log_wgt cmap=cividis \
	-xtitle "Maximum Mass $$ M_\mathrm{max} [M_\odot] $$" \
	-ytitle "Mass Distribution, $$ \mathcal{P}(M=1.6 M_\odot $$)" \
	-ttext 0.5 0.96 "$(model): NS-WD" \
	-save figs/scat_sn-mmax_$(model).png -show &
	o2graph -read tmp/main2_$(model).o2 markov_chain_0 -function 1 w -set colbar 1 \
	-scatter M_max SN_LM w log_wgt cmap=cividis \
	-xtitle "Maximum Mass $$ M_\mathrm{max} [M_\odot] $$" \
	-ytitle "Mass Distribution, $$ \mathcal{P}(M=1.7 M_\odot $$)" \
	-ttext 0.5 0.96 "$(model): LMXB" \
	-save figs/scat_sn-mmax_$(model).png -show &

splot-mmax-sn:
	o2graph -read tmp/main2_$(model).o2 markov_chain_0 -subplots 3 1 sharex=True \
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
	-selax 0 -read tmp/dist_NL.o2 sn_ns -den-plot prob_ns cmap=pink_r \
	-read tmp/dist_NL.o2 c68_ns -plot lw=1,color=purple \
	-read tmp/dist_NL.o2 c95_ns -plot lw=1,color=orange \
	-ttext 0.5 1.05 "DNS" \
	-ttext 0.2 0.92 "NL" \
	-ytitle "Probability Density, $$ \mathcal{P}(M) $$" \
	-selax 1 -read tmp/dist_NL.o2 sn_wd -den-plot prob_wd cmap=pink_r \
	-read tmp/dist_NL.o2 c68_wd -plot lw=1,color=purple \
	-read tmp/dist_NL.o2 c95_wd -plot lw=1,color=orange \
	-ttext 0.5 1.05 "NS-WD" \
	-selax 2 -read tmp/dist_NL.o2 sn_lm -den-plot prob_lm cmap=pink_r \
	-read tmp/dist_NL.o2 c68_lm -plot lw=1,color=purple \
	-read tmp/dist_NL.o2 c95_lm -plot lw=1,color=orange \
	-ttext 0.5 1.05 "LMXB" \
	-selax 3 -read tmp/dist_NP.o2 sn_ns -den-plot prob_ns cmap=pink_r \
	-read tmp/dist_NP.o2 c68_ns -plot lw=1,color=purple \
	-read tmp/dist_NP.o2 c95_ns -plot lw=1,color=orange \
	-ttext 0.2 0.92 "NP" \
	-selax 4 -read tmp/dist_NP.o2 sn_wd -den-plot prob_wd cmap=pink_r \
	-read tmp/dist_NP.o2 c68_wd -plot lw=1,color=purple \
	-read tmp/dist_NP.o2 c95_wd -plot lw=1,color=orange \
	-xtitle "Mass, $$ M $$ [$$ M_\odot $$]" \
	-selax 5 -read tmp/dist_NP.o2 sn_lm -den-plot prob_lm cmap=pink_r \
	-read tmp/dist_NP.o2 c68_lm -plot lw=1,color=purple \
	-read tmp/dist_NP.o2 c95_lm -plot lw=1,color=orange \
	-addcbar 0.94 0.12 0.02 0.83 -selax \
	-save figs/sn_s1.png -show &

splot-dist2:
	o2graph -subplots 2 3 sharex=True,sharey=True \
	-subadj hspace=0.0,left=0.05,right=0.93,top=0.95,bottom=0.12 \
	-selax 0 -read tmp/dist_ML.o2 sn_ns -den-plot prob_ns cmap=pink_r \
	-read tmp/dist_ML.o2 c68_ns -plot lw=1,color=purple \
	-read tmp/dist_ML.o2 c95_ns -plot lw=1,color=orange \
	-ttext 0.5 1.05 "DNS" \
	-ttext 0.2 0.92 "ML" \
	-ytitle "Probability Density, $$ \mathcal{P}(M) $$" \
	-selax 1 -read tmp/dist_ML.o2 sn_wd -den-plot prob_wd cmap=pink_r \
	-read tmp/dist_ML.o2 c68_wd -plot lw=1,color=purple \
	-read tmp/dist_ML.o2 c95_wd -plot lw=1,color=orange \
	-ttext 0.5 1.05 "NS-WD" \
	-selax 2 -read tmp/dist_ML.o2 sn_lm -den-plot prob_lm cmap=pink_r \
	-read tmp/dist_ML.o2 c68_lm -plot lw=1,color=purple \
	-read tmp/dist_ML.o2 c95_lm -plot lw=1,color=orange \
	-ttext 0.5 1.05 "LMXB" \
	-selax 3 -read tmp/dist_MP.o2 sn_ns -den-plot prob_ns cmap=pink_r \
	-read tmp/dist_MP.o2 c68_ns -plot lw=1,color=purple \
	-read tmp/dist_MP.o2 c95_ns -plot lw=1,color=orange \
	-ttext 0.2 0.92 "MP" \
	-selax 4 -read tmp/dist_MP.o2 sn_wd -den-plot prob_wd cmap=pink_r \
	-read tmp/dist_MP.o2 c68_wd -plot lw=1,color=purple \
	-read tmp/dist_MP.o2 c95_wd -plot lw=1,color=orange \
	-xtitle "Mass, $$ M $$ [$$ M_\odot $$]" \
	-selax 5 -read tmp/dist_MP.o2 sn_lm -den-plot prob_lm cmap=pink_r \
	-read tmp/dist_MP.o2 c68_lm -plot lw=1,color=purple \
	-read tmp/dist_MP.o2 c95_lm -plot lw=1,color=orange \
	-addcbar 0.94 0.12 0.02 0.83 -selax \
	-save figs/sn_s2.png -show &

splot-dist:
	o2graph -subplots 4 3 sharex=True,sharey=True \
	-subadj hspace=0.0,left=0.06,right=0.92,top=0.98,bottom=0.05 \
	-selax 0 -read tmp/dist_NL.o2 sn_ns -den-plot prob_ns cmap=pink_r \
	-read tmp/dist_NL.o2 c68_ns -plot lw=1,color=purple \
	-read tmp/dist_NL.o2 c95_ns -plot lw=1,color=orange \
	-ttext 0.2 0.92 "NL" \
	-ttext 0.5 1.05 "DNS" \
	-selax 1 -read tmp/dist_NL.o2 sn_wd -den-plot prob_wd cmap=pink_r \
	-read tmp/dist_NL.o2 c68_wd -plot lw=1,color=purple \
	-read tmp/dist_NL.o2 c95_wd -plot lw=1,color=orange \
	-ttext 0.5 1.05 "NS-WD" \
	-selax 2 -read tmp/dist_NL.o2 sn_lm -den-plot prob_lm cmap=pink_r \
	-read tmp/dist_NL.o2 c68_lm -plot lw=1,color=purple \
	-read tmp/dist_NL.o2 c95_lm -plot lw=1,color=orange \
	-ttext 0.5 1.05 "LMXB" \
	-selax 3 -read tmp/dist_NP.o2 sn_ns -den-plot prob_ns cmap=pink_r \
	-read tmp/dist_NP.o2 c68_ns -plot lw=1,color=purple \
	-read tmp/dist_NP.o2 c95_ns -plot lw=1,color=orange \
	-ttext 0.2 0.92 "NP" \
	-selax 4 -read tmp/dist_NP.o2 sn_wd -den-plot prob_wd cmap=pink_r \
	-read tmp/dist_NP.o2 c68_wd -plot lw=1,color=purple \
	-read tmp/dist_NP.o2 c95_wd -plot lw=1,color=orange \
	-xtitle "Mass, $$ M $$ [$$ M_\odot $$]" \
	-selax 5 -read tmp/dist_NP.o2 sn_lm -den-plot prob_lm cmap=pink_r \
	-read tmp/dist_NP.o2 c68_lm -plot lw=1,color=purple \
	-read tmp/dist_NP.o2 c95_lm -plot lw=1,color=orange \
	-selax 6 -read tmp/dist_ML.o2 sn_ns -den-plot prob_ns cmap=pink_r \
	-read tmp/dist_ML.o2 c68_ns -plot lw=1,color=purple \
	-read tmp/dist_ML.o2 c95_ns -plot lw=1,color=orange \
	-ttext 0.2 0.92 "ML" \
	-ytitle ".$$ \qquad \qquad \qquad \qquad \qquad \qquad \quad $$ Probability Density, $$ \mathcal{P}(M) $$" \
	-selax 7 -read tmp/dist_ML.o2 sn_wd -den-plot prob_wd cmap=pink_r \
	-read tmp/dist_ML.o2 c68_wd -plot lw=1,color=purple \
	-read tmp/dist_ML.o2 c95_wd -plot lw=1,color=orange \
	-selax 8 -read tmp/dist_ML.o2 sn_lm -den-plot prob_lm cmap=pink_r \
	-read tmp/dist_ML.o2 c68_lm -plot lw=1,color=purple \
	-read tmp/dist_ML.o2 c95_lm -plot lw=1,color=orange \
	-selax 9 -read tmp/dist_MP.o2 sn_ns -den-plot prob_ns cmap=pink_r \
	-read tmp/dist_MP.o2 c68_ns -plot lw=1,color=purple \
	-read tmp/dist_MP.o2 c95_ns -plot lw=1,color=orange \
	-ttext 0.2 0.92 "MP" \
	-selax 10 -read tmp/dist_MP.o2 sn_wd -den-plot prob_wd cmap=pink_r \
	-read tmp/dist_MP.o2 c68_wd -plot lw=1,color=purple \
	-read tmp/dist_MP.o2 c95_wd -plot lw=1,color=orange \
	-xtitle "Mass, $$ M $$ [$$ M_\odot $$]" \
	-selax 11 -read tmp/dist_MP.o2 sn_lm -den-plot prob_lm cmap=pink_r \
	-read tmp/dist_MP.o2 c68_lm -plot lw=1,color=purple \
	-read tmp/dist_MP.o2 c95_lm -plot lw=1,color=orange \
	-addcbar 0.93 0.05 0.02 0.93 -selax \
	-save figs/sn_s.png -show &

plot-dist:
	o2graph -read tmp/dist_$(model).o2 sn_ns -set colbar 1 \
	-den-plot prob_ns cmap=pink_r \
	-read tmp/dist_$(model).o2 c68_ns -plot lw=1,color=purple \
	-read tmp/dist_$(model).o2 c95_ns -plot lw=1,color=orange \
	-xtitle "Mass, $$ M $$ [$$ M_\odot $$]" \
	-ytitle "Mass Distribution, $$ \mathcal{P} $$" \
	-ttext 0.5 0.96 "$(model): DNS" \
	-save figs/sn_ns_$(model).png -show & 
	o2graph -read tmp/dist_$(model).o2 sn_wd -set colbar 1 \
	-den-plot prob_wd cmap=pink_r \
	-read tmp/dist_$(model).o2 c68_wd -plot lw=1,color=purple \
	-read tmp/dist_$(model).o2 c95_wd -plot lw=1,color=orange \
	-xtitle "Mass, $$ M $$ [$$ M_\odot $$]" \
	-ytitle "Mass Distribution, $$ \mathcal{P} $$" \
	-ttext 0.5 0.96 "$(model): NS-WD" \
	-save figs/sn_wd_$(model).png -show & 
	o2graph -read tmp/dist_$(model).o2 sn_lm -set colbar 1 \
	-den-plot prob_lm cmap=pink_r \
	-read tmp/dist_$(model).o2 c68_lm -plot lw=1,color=purple \
	-read tmp/dist_$(model).o2 c95_lm -plot lw=1,color=orange \
	-xtitle "Mass, $$ M $$ [$$ M_\odot $$]" \
	-ytitle "Mass Distribution, $$ \mathcal{P} $$" \
	-ttext 0.5 0.96 "$(model): LMXB" \
	-save figs/sn_lm_$(model).png -show & 

file-eos:
	acol -read tmp/main_$(model).o2 markov_chain_0 \
	$(foreach i,$(INDEX),-function P_$i*197.32698 Phc_$i) \
	-set obj_name ed_pr -select ":Phc_*" -internal tmp/eos_$(model).o2 \
	-read tmp/eos_$(model).o2 ed_pr -ser-hist-t3d \
	"grid:0.3*197.32698,10.0*197.32698,((10.0-0.3)*197.32698)/99.0" x ed \
	"grid:1.0,2.0e3,(2.0e3)^(1.0/99.0),log" auto pr \
	"Phc_*" prob -type -list \
	-set obj_name ed_pr_3d -internal tmp/eos_$(model).o2 \
	-to-hist-2d prob \
	-contours frac 0.68 tmp/eos_$(model).o2 c68 \
	-contours frac 0.95 tmp/eos_$(model).o2 c95

plot-eos:
	o2graph -read tmp/eos_$(model).o2 ed_pr_3d -set logy 1 -set colbar 1 \
	-den-plot prob pcm=True,cmap=pink_r \
	-read tmp/eos_$(model).o2 c68 -plot lw=1,color=purple \
	-read tmp/eos_$(model).o2 c95 -plot lw=1,color=orange \
	-xtitle "Energy Density, $$ \epsilon~[\mathrm{MeV}/\mathrm{fm}^3] $$" \
	-ytitle "Pressure, $$ P~[\mathrm{MeV}/\mathrm{fm}^3] $$" \
	-ttext 0.5 0.96 "$(model): EoS" \
	-save figs/eos_$(model).png -show &

file-mvsr:
	acol -read tmp/main_$(model).o2 markov_chain_0 -delete-col R_max \
	-set obj_name gm_r -select ":R_*" -internal tmp/mvsr_$(model).o2 \
	-read tmp/mvsr_$(model).o2 gm_r -ser-hist-t3d \
	"grid:0.2,3.0,(3.0-0.2)/99.0" y gm \
	"grid:9.0,15.0,(15.0-9.0)/99.0" auto r \
	"R_*" prob -type -list \
	-set obj_name gm_r_3d -internal tmp/mvsr_$(model).o2 \
	-to-hist-2d prob \
	-contours frac 0.68 tmp/mvsr_$(model).o2 c68 \
	-contours frac 0.95 tmp/mvsr_$(model).o2 c95

plot-mvsr:
	o2graph -read tmp/mvsr_$(model).o2 gm_r_3d -set colbar 1 \
	-den-plot prob cmap=pink_r \
	-read tmp/mvsr_$(model).o2 c68 -plot lw=1,color=purple \
	-read tmp/mvsr_$(model).o2 c95 -plot lw=1,color=orange \
	-xtitle "Radius, $$ R~[\mathrm{km}] $$" \
	-ytitle "Mass, $$ M~[\mathrm{M}_\odot] $$" \
	-ttext 0.5 0.96 "$(model): $(title)" \
	-save figs/mvsr_$(model).png -show &

file-eoa:
	acol -read tmp/main_$(model).o2 markov_chain_0 -select ":EoA_*" \
	-set obj_name nb_ea -internal tmp/eoa_$(model).o2 \
	-read tmp/eoa_$(model).o2 nb_ea -ser-hist-t3d \
	"grid:0.04,1.24,(1.24-0.04)/99.0" x nb \
	"grid:0.0,2.5,(2.5-0.0)/99.0" auto ea \
	"EoA_*" prob -type -list \
	-set obj_name nb_ea_3d -internal tmp/eoa_$(model).o2 \
	-to-hist-2d prob \
	-contours frac 0.68 tmp/eoa_$(model).o2 c68 \
	-contours frac 0.95 tmp/eoa_$(model).o2 c95 \
	-contours frac 0.99 tmp/eoa_$(model).o2 c99

plot-eoa:
	o2graph -read tmp/eoa_$(model).o2 nb_ea_3d -set colbar 1 \
	-den-plot prob cmap=pink_r \
	-read tmp/eoa_$(model).o2 c68 -plot lw=1,color=purple \
	-read tmp/eoa_$(model).o2 c95 -plot lw=1,color=orange \
	-xtitle "Baryon Density, $$ n_b~[\mathrm{fm}^{-3}] $$" \
	-ytitle "Energy/Nucleon, $$ E/A~[\mathrm{MeV}] $$" \
	-ttext 0.5 0.96 "$(model): $(title)" \
	-save figs/eoa_$(model).png -show &

file-pnb:
	acol -read tmp/main_$(model).o2 markov_chain_0 \
	$(foreach i,$(INDEX),-function Pnb_$i*197.32698 Pnb_hc_$i) \
	-select ":Pnb_hc_*" -set obj_name nb_pr -internal tmp/pnb_$(model).o2 \
	-read tmp/pnb_$(model).o2 nb_pr -ser-hist-t3d \
	"grid:0.04,1.24,(1.24-0.04)/99.0" x nb \
	"grid:1.0,2.0e3,(2.0e3)^(1/99),log" auto pr \
	"Pnb_hc_*" prob -type -list \
	-set obj_name nb_pr_3d -internal tmp/pnb_$(model).o2 \
	-to-hist-2d prob \
	-contours frac 0.68 tmp/pnb_$(model).o2 c68 \
	-contours frac 0.95 tmp/pnb_$(model).o2 c95 \
	-contours frac 0.99 tmp/pnb_$(model).o2 c99

plot-pnb:
	o2graph -read tmp/pnb_$(model).o2 nb_pr_3d \
	-set colbar 1 -set logy 1 \
	-den-plot prob pcm=True,cmap=pink_r \
	-read tmp/pnb_$(model).o2 c68 -plot lw=1,color=purple \
	-read tmp/pnb_$(model).o2 c95 -plot lw=1,color=orange \
	-xtitle "Baryon Density, $$ n_b~[\mathrm{fm}^{-3}] $$" \
	-ytitle "Pressure, $$ P~[\mathrm{MeV/fm}^3] $$" \
	-ttext 0.5 0.96 "$(model): $(title)" \
	-save figs/pnb_$(model).png -show &

plot-pnb2:
	o2graph -plotv "grid:0.04,1.24,(1.24-0.04)/99.0" \
	"hdf5:tmp/pnb_ML.o2:nb_pr:0-99:Pnb_hc_*" \
	-xtitle "Baryon Density, $ n_b~[\mathrm{fm}^{-3}] $" \
	-ytitle "Pressure, $ P~[\mathrm{MeV/fm}^3] $" -show &

plot-trc-wgt:
	o2graph -read tmp/main_$(model).o2 markov_chain_0 \
	-xtitle "Steps" -ytitle "Average log_wgt" \
	-ttext 0.5 0.96 "$(model): $(title)" \
	-average-rows log_wgt 512 -plot1 log_wgt \
	-save figs/trc_wgt_$(model).png -show &

plot-trc-mmax:
	o2graph -read tmp/main_$(model).o2 markov_chain_0 \
	-xtitle "Steps" -ytitle "Average M_max" \
	-ttext 0.5 0.96 "$(model): $(title)" \
	-average-rows M_max 512 -plot1 M_max \
	-save figs/trc_mmax_$(model).png -show &

plot-trc-a:
	o2graph -read tmp/main_$(model).o2 markov_chain_0 \
	-xtitle "Steps" -ytitle "Average $$ a $$ [MeV]" \
	-ttext 0.5 0.96 "$(model): $(title)" \
	-average-rows a 512 -plot1 a \
	-save figs/trc_a_$(model).png -show &

plot-trc-alpha:
	o2graph -read tmp/main_$(model).o2 markov_chain_0 \
	-xtitle "Steps" -ytitle "Average $$ \alpha $$ [MeV]" \
	-ttext 0.5 0.96 "$(model): $(title)" \
	-average-rows alpha 512 -plot1 alpha \
	-save figs/trc_alpha_$(model).png -show &

plot-trc-r43:
	o2graph -read tmp/main_$(model).o2 markov_chain_0 \
	-xtitle "Steps" -ytitle "Average $$ R_{1.4M_\odot} $$ [km]" \
	-ttext 0.5 0.96 "$(model): $(title)" \
	-average-rows R_43 512 -plot1 R_43 \
	-save figs/trc_r43_$(model).png -show &

plot-trc-dpdm:
	o2graph -read tmp/main_$(model).o2 markov_chain_0 \
	-function "log(dpdM)" log_dpdm \
	-xtitle "Steps" -ytitle "Average log(dp/dM)" \
	-ttext 0.5 0.96 "$(model): $(title)" \
	-average-rows log_dpdm 512 -plot1 log_dpdm \
	-save figs/trc_dpdm_$(model).png -show &

plot-scat-mmax-wgt:
	o2graph -read tmp/main_$(model).o2 markov_chain_0 \
	-function 1 w -scatter M_max log_wgt w log_wgt \
	-xtitle "M_max [$$ M_\odot $$]" -ytitle "log_wgt" \
	-ttext 0.5 0.96 "$(model): $(title)" \
	-save figs/scat_mmax-wgt_$(model).png -show &

plot-scat-dpdm-mmax:
	o2graph -read tmp/main_$(model).o2 markov_chain_0 \
	-function "log(dpdM)" log_dpdm -set colbar 1 \
	-function 1 w -scatter M_max log_dpdm w log_wgt \
	-xtitle "M_max [$$ M_\odot $$]" -ytitle "log(dpdM)" \
	-ttext 0.5 0.96 "$(model): $(title)" \
	-save figs/scat_dpdm-mmax_$(model).png -show &

plot-scat-dpdm-poly:
	o2graph -read tmp/main_$(model).o2 markov_chain_0 \
	-function "log(dpdM)" log_dpdm -set colbar 1 \
	-function 1 w -scatter exp3 log_dpdm w log_wgt \
	-xtitle "exp3" -ytitle "log(dp/dM)" \
	-ttext 0.5 0.96 "$(model): $(title)" \
	-save figs/scat_dpdm-poly_$(model).png -show &

plot-scat-dpdm-line:
	o2graph -read tmp/main_$(model).o2 markov_chain_0 \
	-function "log(dpdM)" log_dpdm -set colbar 1 \
	-function 1 w -scatter csq3 log_dpdm w log_wgt \
	-xtitle "csq3" -ytitle "log(dp/dM)" \
	-ttext 0.5 0.96 "$(model): $(title)" \
	-save figs/scat_dpdm-line_$(model).png -show &

plot-scat-dpdm-r43:
	o2graph -read tmp/main_$(model).o2 markov_chain_0 \
	-function "log(dpdM)" log_dpdm -set colbar 1 \
	-function 1 w -scatter R_43 log_dpdm w log_wgt \
	-xtitle "R_43 [km]" -ytitle "log(dpdM)" \
	-ttext 0.5 0.96 "$(model): $(title)" \
	-save figs/scat_dpdm-r43_$(model).png -show &

plot-acor-a:
	o2graph $(foreach i, $(THREAD),$(foreach j,$(WALKER), \
	-read tmp/main_$(model).o2 markov_chain_0 \
	-select-rows "thread==$i && walker==$j" \
	-autocorr store,default a -plot 0)) \
	-xtitle "Autocorr. Length" -ytitle "Autocorr. Coeff." \
	-ttext 0.5 0.96 "$(model): $(title)" -save figs/a-acor_$(model).png -show &

plot-kde-r43:
	o2graph -read mp334 markov_chain_0 -kde-plot R_43 -read mp356 markov_chain_0 -kde-plot R_43 \
	-xtitle "R_43" -ytitle "Density" -ttext 0.5 0.97 "MP: 33-34 vs 35-36" -save r43_mp.png -show &
	o2graph -read nl334 markov_chain_0 -kde-plot R_43 -read nl356 markov_chain_0 -kde-plot R_43 \
	-xtitle "R_43" -ytitle "Density" -ttext 0.5 0.97 "NL: 33-34 vs 35-36" -save r43_nl.png -show &
	o2graph -read np334 markov_chain_0 -kde-plot R_43 -read np356 markov_chain_0 -kde-plot R_43 \
	-xtitle "R_43" -ytitle "Density" -ttext 0.5 0.97 "NP: 33-34 vs 35-36" -save r43_np.png -show &

plot-kde-mmax:
	o2graph -read mp334 markov_chain_0 -kde-plot M_max -read mp356 markov_chain_0 -kde-plot M_max \
	-xtitle "M_max" -ytitle "Density" -ttext 0.5 0.97 "MP: 33-34 vs 35-36" -save mmax_mp.png -show &
	o2graph -read nl334 markov_chain_0 -kde-plot M_max -read nl356 markov_chain_0 -kde-plot M_max \
	-xtitle "M_max" -ytitle "Density" -ttext 0.5 0.97 "NL: 33-34 vs 35-36" -save mmax_nl.png -show &
	o2graph -read np334 markov_chain_0 -kde-plot M_max -read np356 markov_chain_0 -kde-plot M_max \
	-xtitle "M_max" -ytitle "Density" -ttext 0.5 0.97 "NP: 33-34 vs 35-36" -save mmax_np.png -show &

plot-kde-a:
	o2graph -read mp334 markov_chain_0 -kde-plot a -read mp356 markov_chain_0 -kde-plot a \
	-xtitle "a" -ytitle "Density" -ttext 0.5 0.97 "MP: 33-34 vs 35-36" -save a_mp.png -show &
	o2graph -read nl334 markov_chain_0 -kde-plot a -read nl356 markov_chain_0 -kde-plot a \
	-xtitle "a" -ytitle "Density" -ttext 0.5 0.97 "NL: 33-34 vs 35-36" -save a_nl.png -show &
	o2graph -read np334 markov_chain_0 -kde-plot a -read np356 markov_chain_0 -kde-plot a \
	-xtitle "a" -ytitle "Density" -ttext 0.5 0.97 "NP: 33-34 vs 35-36" -save a_np.png -show &

plot-kde-q:
	o2graph -read mp334 markov_chain_0 -kde-plot q -read mp356 markov_chain_0 -kde-plot q \
	-xtitle "q" -ytitle "Density" -ttext 0.5 0.97 "MP: 33-34 vs 35-36" -save q_mp.png -show &
	o2graph -read nl334 markov_chain_0 -kde-plot q -read nl356 markov_chain_0 -kde-plot q \
	-xtitle "q" -ytitle "Density" -ttext 0.5 0.97 "NL: 33-34 vs 35-36" -save q_nl.png -show &
	o2graph -read np334 markov_chain_0 -kde-plot q -read np356 markov_chain_0 -kde-plot q \
	-xtitle "q" -ytitle "Density" -ttext 0.5 0.97 "NP: 33-34 vs 35-36" -save q_np.png -show &

plot-kde-m1gw19:
	o2graph -read mp334 markov_chain_0 -kde-plot m1_gw19 -read mp356 markov_chain_0 -kde-plot m1_gw19 \
	-xtitle "m1_gw19" -ytitle "Density" -ttext 0.5 0.97 "MP: 33-34 vs 35-36" -save m1gw19_mp.png -show &
	o2graph -read nl334 markov_chain_0 -kde-plot m1_gw19 -read nl356 markov_chain_0 -kde-plot m1_gw19 \
	-xtitle "m1_gw19" -ytitle "Density" -ttext 0.5 0.97 "NL: 33-34 vs 35-36" -save m1gw19_nl.png -show &
	o2graph -read np334 markov_chain_0 -kde-plot m1_gw19 -read np356 markov_chain_0 -kde-plot m1_gw19 \
	-xtitle "m1_gw19" -ytitle "Density" -ttext 0.5 0.97 "NP: 33-34 vs 35-36" -save m1gw19_np.png -show &

plot-kde-wgt:
	o2graph -read mp334 markov_chain_0 -kde-plot log_wgt -read mp356 markov_chain_0 -kde-plot log_wgt \
	-xtitle "log_wgt" -ytitle "Density" -ttext 0.5 0.97 "MP: 33-34 vs 35-36" -save wgt_mp.png -show &
	o2graph -read nl334 markov_chain_0 -kde-plot log_wgt -read nl356 markov_chain_0 -kde-plot log_wgt \
	-xtitle "log_wgt" -ytitle "Density" -ttext 0.5 0.97 "NL: 33-34 vs 35-36" -save wgt_nl.png -show &
	o2graph -read np334 markov_chain_0 -kde-plot log_wgt -read np356 markov_chain_0 -kde-plot log_wgt \
	-xtitle "log_wgt" -ytitle "Density" -ttext 0.5 0.97 "NP: 33-34 vs 35-36" -save wgt_np.png -show &

plot-kde-mean:
	o2graph -read mp334 markov_chain_0 -kde-plot mean_NS -read mp356 markov_chain_0 -kde-plot mean_NS \
	-xtitle "mean_NS" -ytitle "Density" -ttext 0.5 0.97 "MP: 33-34 vs 35-36" -save mean_mp.png -show &
	o2graph -read nl334 markov_chain_0 -kde-plot mean_NS -read nl356 markov_chain_0 -kde-plot mean_NS \
	-xtitle "mean_NS" -ytitle "Density" -ttext 0.5 0.97 "NL: 33-34 vs 35-36" -save mean_nl.png -show &
	o2graph -read np334 markov_chain_0 -kde-plot mean_NS -read np356 markov_chain_0 -kde-plot mean_NS \
	-xtitle "mean_NS" -ytitle "Density" -ttext 0.5 0.97 "NP: 33-34 vs 35-36" -save mean_np.png -show &

plot-kde-dpdm:
	o2graph -read mp334 markov_chain_0 -function "log(dpdM)" log_dpdm -kde-plot log_dpdm -read mp356 markov_chain_0 -function "log(dpdM)" log_dpdm -kde-plot log_dpdm \
	-xtitle "log(dp/dM)" -ytitle "Density" -ttext 0.5 0.97 "MP: 33-34 vs 35-36" -save dpdm_mp.png -show &
	o2graph -read nl334 markov_chain_0 -function "log(dpdM)" log_dpdm -kde-plot log_dpdm -read nl356 markov_chain_0 -function "log(dpdM)" log_dpdm -kde-plot log_dpdm \
	-xtitle "log(dp/dM)" -ytitle "Density" -ttext 0.5 0.97 "NL: 33-34 vs 35-36" -save dpdm_nl.png -show &
	o2graph -read np334 markov_chain_0 -function "log(dpdM)" log_dpdm -kde-plot log_dpdm -read np356 markov_chain_0 -function "log(dpdM)" log_dpdm -kde-plot log_dpdm \
	-xtitle "log(dp/dM)" -ytitle "Density" -ttext 0.5 0.97 "NP: 33-34 vs 35-36" -save dpdm_np.png -show &

plot-kde-csq:
	o2graph -read ml334 markov_chain_0 -function "csq1-csq2" dcsq -kde-plot dcsq -read ml356 markov_chain_0 -function "csq1-csq2" dcsq -kde-plot dcsq \
	-xtitle "csq1-csq2" -ytitle "Density" -ttext 0.5 0.97 "ML: 33-34 vs 35-36" -save csq_ml.png -show &
	o2graph -read nl334 markov_chain_0 -function "csq1-csq2" dcsq -kde-plot dcsq -read nl356 markov_chain_0 -function "csq1-csq2" dcsq -kde-plot dcsq \
	-xtitle "csq1-csq2" -ytitle "Density" -ttext 0.5 0.97 "NL: 33-34 vs 35-36" -save csq_nl.png -show &

plot-r43-lambda43:
	o2graph -read np_43_tmp markov_chain_0 -index \
	i=1; while [[$$i -le 8034]]; do \
		-select-rows "N==$$i" \
		-scatter R_43 lambda_43; \
		i=$$((i+1)); \
	done
	

file-eos-all:
	make model="ML" file-eos && \
	make model="MP" file-eos && \
	make model="NL" file-eos && \
	make model="NP" file-eos

file-eoa-all:
	make model="ML" file-eoa && \
	make model="MP" file-eoa && \
	make model="NL" file-eoa && \
	make model="NP" file-eoa

file-pnb-all:
	make model="ML" file-pnb && \
	make model="MP" file-pnb && \
	make model="NL" file-pnb && \
	make model="NP" file-pnb

file-mvsr-all:
	make model="ML" file-mvsr && \
	make model="MP" file-mvsr && \
	make model="NL" file-mvsr && \
	make model="NP" file-mvsr

plot-eos-all:
	make model="ML" plot-eos && \
	make model="MP" plot-eos && \
	make model="NL" plot-eos && \
	make model="NP" plot-eos

plot-eoa-all:
	make model="ML" plot-eoa && \
	make model="MP" plot-eoa && \
	make model="NL" plot-eoa && \
	make model="NP" plot-eoa

plot-pnb-all:
	make model="ML" plot-pnb && \
	make model="MP" plot-pnb && \
	make model="NL" plot-pnb && \
	make model="NP" plot-pnb

plot-mvsr-all:
	make model="ML" plot-mvsr && \
	make model="MP" plot-mvsr && \
	make model="NL" plot-mvsr && \
	make model="NP" plot-mvsr

file-all:
	make model="ML" file-esmr && make model="MP" file-esmr && make model="NL" file-esmr && make model="NP" file-esmr && \
	make model="ML" file-main2 && make model="MP" file-main2 && make model="NL" file-main2 && make model="NP" file-main2 && \
	make model="ML" file-dist && make model="MP" file-dist && make model="NL" file-dist && make model="NP" file-dist

file-dist-all:
	make model="ML" file-dist && make model="MP" file-dist && make model="NL" file-dist && make model="NP" file-dist

plot-dist-all:
	make model="ML" plot-dist && make model="MP" plot-dist && make model="NL" plot-dist && make model="NP" plot-dist

plot-all:
	make plot-eos && \
	make plot-mvsr && \
	make plot-dist && \
	make plot-trc && \
	make plot-trc-mmax && \
	make plot-scat-mmax && \
	make plot-scat-dpdm-mmax && \
	make plot-scat-dpdm-r43

clean:
	rm -f figs/*.png tmp/*.o2 