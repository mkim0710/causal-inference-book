 /***************************************************************
PROGRAM 15.1
Estimating the average causal effect within levels of confounders
under the assumption of effect-measure modification by smoking
intensity ONLY 
Data from NHEFS
***************************************************************/

clear
use "P:\CIBook\NHEFSdata\nhefs.dta"

gen byte cens = (wt82 == .)

/* provisionally ignore subjects with missing values*/
*Analysis restricted to N=1566 with no missing values*/
foreach var of varlist wt82 education{
drop if `var'==.
}

* regression on covariates, allowing for some effect modfication
regress wt82_71 c.qsmk##c.smokeintensity c.smokeintensity##c.smokeintensity sex race c.age##c.age ///
ib(last).education c.smokeyrs##c.smokeyrs ib(last).exercise ib(last).active c.wt71##c.wt71 



/* regression on covaraites, not allowing for effect modification*/
regress wt82_71 qsmk c.smokeintensity##c.smokeintensity sex race c.age##c.age ///
ib(last).education c.smokeyrs##c.smokeyrs ib(last).exercise ib(last).active c.wt71##c.wt71 



/***************************************************************
PROGRAM 15.2
Estimating and plotting the propensity score
Data from NHEFS
***************************************************************/
clear
use "P:\CIBook\NHEFSdata\nhefs.dta"

gen byte cens = (wt82 == .)

/* provisionally ignore subjects with missing values*/
/*Analysis restricted to N=1629 with no missing covariates*/
/*leave missing wt81_72 in the data*/
foreach var of varlist education{
drop if `var'==.
}

logit qsmk sex race c.age##c.age ib(last).education c.smokeintensity##c.smokeintensity ///
c.smokeyrs##c.smokeyrs ib(last).exercise ib(last).active c.wt71##c.wt71 

predict ps, pr

bys qsmk: summarize ps 

/*plotting the estimated propensity score*/

histogram ps, width(0.05) start(0.025) frequency fcolor(none) lcolor(black) lpattern(solid) addlabel addlabopts(mlabcolor(black) mlabposition(12) mlabangle(zero)) ///
ytitle(No. Subjects) ylabel(#4) xtitle(Estimated Propensity Score) xlabel(#15)  ///
by(, title(Estimated Propensity Score Distribution) subtitle(By Quit Smoking Status)) ///
by(, legend(off)) by(qsmk, style(compact)  colfirst) subtitle(, size(small) box bexpand)


/***************************************************************************
PROGRAM 15.3
Stratification and outcome regression using deciles of the propensity score
Data from NHEFS
***************************************************************************/

/*calculation of deciles of ps*/
xtile ps_dec = ps, nq(10)

by ps_dec qsmk, sort: summarize ps

/*stratification on PS deciles, allowing for effect modification*/

by ps_dec: ttest wt82_71, by(qsmk)

/*regression on PS deciles, not allowing for effect modifciation*/
regress wt82_71 qsmk ib(last).ps_dec



/***************************************************************
PROGRAM 15.4
Standardization and outcome regression using the propensity score
Data from NHEFS
***************************************************************/
clear
use "P:\CIBook\NHEFSdata\nhefs.dta"

gen byte cens = (wt82 == .)
/*only discard individuals with missing covariates for PS estimation*/
drop if education==.

logit qsmk sex race c.age##c.age ib(last).education c.smokeintensity##c.smokeintensity ///
c.smokeyrs##c.smokeyrs ib(last).exercise ib(last).active c.wt71##c.wt71 

predict ps, pr
gen meanY_b = .
save nhefs_ps

/*regression on the propensity score, not allowing for effect modification*/
regress wt82_71 qsmk ps

/*standardization by propensity score, agnostic regarding effect modification*/
gen copy1 = 2
expand copy1, generate(interv)
replace interv = -1 if interv == 0
replace interv = 0  if interv == 1
replace copy1 = 1   if interv == 0
expand copy1, generate(interv2)
replace interv = 1  if interv2 ==1
drop interv2 
drop copy1
tab interv
replace wt82_71 = . if interv != -1
replace qsmk = 0 if interv == 0
replace qsmk = 1 if interv == 1
by interv, sort: summarize qsmk


/*estimation in original sample*/
quietly regress wt82_71 qsmk ps
predict predY

by interv, sort: summarize predY

/*store original sample results in a matrix for later comparison*/
gen meanY=.
summarize predY if(interv == -1)
replace meanY = r(mean) if(interv==-1)
summarize meanY if(interv==-1)
matrix input observe = (-1,`r(mean)')
summarize predY if(interv == 0)
replace meanY = r(mean) if(interv==0)
summarize meanY if(interv==0)
matrix observe = (observe \0,`r(mean)')
summarize predY if(interv == 1)
replace meanY = r(mean) if(interv==1)
summarize meanY if(interv==1)
matrix observe = (observe \1,`r(mean)')



/*check that matrix is correct*/
matrix list observe
by interv, sort: summarize predY

program drop bootstdz

/*bootstrap program*/

capture program drop bootstdz
program define bootstdz, rclass
	u nhefs_ps, clear
		preserve
		/*draw bootstrap sample from original observations*/
		bsample 
		/*create copies with each value of qsmk in bootstrap sample*/
		gen copy1 = 2
		expand copy1, generate(interv)
		replace interv = -1 if interv == 0
		replace interv = 0  if interv == 1
		replace copy1 = 1   if interv == 0
		expand copy1, generate(interv2)
		replace interv = 1  if interv2 ==1
		drop interv2 
		drop copy1
		replace wt82_71 = . if interv != -1
		replace qsmk = 0 if interv == 0
		replace qsmk = 1 if interv == 1
		tab interv
		/*run regression*/
		regress wt82_71 qsmk ps
		predict predY_b
		summarize predY_b if(interv == -1)
		replace meanY_b = r(mean) if(interv==-1)
		summarize meanY_b if(interv==-1)
		return scalar boot_obs = r(mean)
		summarize predY_b if(interv == 0)
		replace meanY_b = r(mean) if(interv==0)
		summarize meanY_b if(interv==0)
		return scalar boot_0 = r(mean)
		summarize predY_b if(interv == 1)
		replace meanY_b = r(mean) if(interv==1)
		summarize meanY_b if(interv==1)
		return scalar boot_1 = r(mean)
	restore
end
simulate boot_obs=r(boot_obs), reps(500) seed(1): bootstdz
gen interv = -1
rename boot_obs meanY_b
save boot_obs
simulate boot_0=r(boot_0), reps(500) seed(1): bootstdz
gen interv = 0
rename boot_0 meanY_b
save boot_0
simulate boot_1=r(boot_1), reps(500) seed(1): bootstdz
gen interv = 1
rename boot_1 meanY_b
save boot_1

use boot_obs
append using boot_0
append using boot_1

gen std=.
summarize meanY_b if(interv == -1)
replace std = r(sd) if(interv==-1)
matrix input boot_std = (-1,`r(sd)')
summarize meanY_b if(interv == 0)
replace std = r(sd) if(interv==0)
matrix boot_std= (boot_std \0,`r(sd)')
summarize meanY_b if(interv == 1)
replace std = r(sd) if(interv==1)
matrix boot_std = (boot_std \1,`r(sd)')
matrix list boot_std
matrix list observe

matrix lb = J(3,2,0)
matrix list lb
matrix ub = J(3,2,0)
matrix list ub
forvalues i = 1/3{
	matrix lb[`i',2]	= (observe[`i',2]-1.96*boot_std[`i',2]) 
	matrix lb[`i',1]    = observe[`i',1]
	matrix ub[`i',2]	= (observe[`i',2]+1.96*boot_std[`i',2])
	matrix ub[`i',1]    = observe[`i',1]
}

matrix list lb
matrix list ub
matrix lb = lb[1...,2]
matrix ub = ub[1...,2]

matrix results=observe,lb,ub
matrix list results

