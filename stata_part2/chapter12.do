 /***************************************************************
PROGRAM 12.1
Descriptive statistics from NHEFS data (Table 12.1)
***************************************************************/

clear
use "P:\CIBook\NHEFSdata\nhefs.dta"

gen byte cens = (wt82 == .)
gen byte older = (age > 50 & age < .)
label var older "50+ years"

/* provisionally ignore subjects with missing values*/
*Analysis restricted to N=1566 with no missing values*/
foreach var of varlist school wt82 {
drop if `var'==.
}

/* mean weight change in those with and without smoking cessation*/
label define qsmk 0 "No smoking cessation" 1 "Smoking cessation"
label values qsmk qsmk
label var wt82_71 "weight change in kg"
bysort qsmk: sum wt82_71 if cens == 0
regress wt82_71 qsmk if cens == 0, noheader
by qsmk, sort: egen years = mean(age) if age < . & cens == 0
label var years "Age, years"
by qsmk, sort: egen cigs = mean(smokeintensity) if smokeintensity < . & cens== 0
label var cigs "Cigarettes/day"
by qsmk, sort: egen kg = mean(wt71) if wt71 < . & cens == 0
label var kg "Weight, kg"
by qsmk, sort: egen white = mean(100 * (race==0)) if race < . & cens == 0
label var white "White, %"
by qsmk, sort: egen male = mean(100 * (sex==0)) if sex < . & cens == 0
label var male "Men, %"
by qsmk, sort: egen inactive = mean(100 * (active==2)) if active < . & cens ==0
label var inactive "Inactive daily life"
by qsmk, sort: egen noexer = mean(100 * (exercise == 2)) if exercise < . & cens == 0
label var noexer "Little/no exercise"
by qsmk, sort: egen university = mean(100 * (education == 5)) if education <. & cens == 0
label var university "University education"
by qsmk, sort: egen syrs = mean(smokeyrs) if smokeyrs < . & cens== 0
label var syrs "Years of smoking"
by qsmk, sort: egen fiftyplus = mean(100 * (older==1)) if sex < . & cens == 0
label var fiftyplus "50+ years, %"

/*Table*/
foreach var of varlist years male white university kg cigs noexer inactive syrs fiftyplus {
tabdisp qsmk, cell(`var') format(%3.1f)
}



/***************************************************************
PROGRAM 12.2
Estimating IP weights
Data from NHEFS
***************************************************************/

* estimation of ip weights via a logistic model 
logit qsmk sex race c.age##c.age ib(last).education c.smokeintensity##c.smokeintensity ///
c.smokeyrs##c.smokeyrs ib(last).exercise ib(last).active c.wt71##c.wt71 
predict p_qsmk, pr

gen w=.
replace w=1/p_qsmk if qsmk==1
replace w=1/(1-p_qsmk) if qsmk==0
summarize w

regress wt82_71 qsmk [pweight=w], cluster(seqn) 


/* no association between sex and qsmk in pseudo-population*/
tab sex qsmk [aw=w]

/* "check" for positivity*/ 
tab  age qsmk if race==0 & sex==1 & wt82!=.



/***************************************************************
PROGRAM 12.3
Estimating stabilized IP weights
Data from NHEFS
***************************************************************/

/* estimation of denominator of ip weights*/ 
logit qsmk sex race c.age##c.age ib(last).education c.smokeintensity##c.smokeintensity ///
c.smokeyrs##c.smokeyrs ib(last).exercise ib(last).active c.wt71##c.wt71 

predict pd_qsmk, pr

/* estimation of numerator of ip weights*/ 
logit qsmk 
predict pn_qsmk, pr

gen sw_a=.
replace sw_a=pn_qsmk/pd_qsmk if qsmk==1
replace sw_a=(1-pn_qsmk)/(1-pd_qsmk) if qsmk==0
summarize sw_a

regress wt82_71 qsmk [pweight=sw_a], cluster(seqn) 

/* no association between sex and qsmk in pseudo-population*/
tab sex qsmk [aw=sw_a]



/***************************************************************
PROGRAM 12.4
Estimating the parameters of a marginal structural mean model
with a continuous treatment Data from NHEFS
***************************************************************/

clear
use "P:\CIBook\NHEFSdata\nhefs.dta"


/* provisionally ignore subjects with missing values*/ 
foreach var of varlist qsmk sex race age school smokeintensity smokeyrs exercise active wt71 wt82 {
drop if `var'==.
}

/* Analysis restricted to subjects reporting <=25 cig/day at baseline*/
keep if  smokeintensity <=25


/* estimation of denominator of ip weights*/

regress smkintensity82_71 sex race c.age##c.age ib(last).education c.smokeintensity##c.smokeintensity ///
c.smokeyrs##c.smokeyrs ib(last).exercise ib(last).active c.wt71##c.wt71

quietly predict p_den
gen dens_den = normalden(smkintensity82_71, p_den, e(rmse))
*NOTE: The regress command in STATA saves the root mean squared error for the immediate regression as e(rmse), thus there is no need to calculate it again. 

/* estimation of numerator of ip weights*/

quietly regress smkintensity82_71
quietly predict p_num

gen dens_num = normalden( smkintensity82_71, p_num, e(rmse))

/* estimation of Stabilized weights*/
gen sw_a=dens_num/dens_den
summarize sw_a
regress wt82_71  c.smkintensity82_71##c. smkintensity82_71 [pweight=sw_a], cluster(seqn)



/***************************************************************
PROGRAM 12.5
Estimating the parameters of a marginal structural logistic model
Data from NHEFS
***************************************************************/
clear
use "P:\CIBook\NHEFSdata\nhefs.dta"


* provisionally ignore subjects with missing values 
foreach var of varlist qsmk sex race age school smokeintensity smokeyrs exercise active wt71 wt82 {
drop if `var'==.
}

* First, estimation of stabilized weights sw_a (same as in PROGRAM 12.3)
logit qsmk sex race c.age##c.age ib(last).education c.smokeintensity##c.smokeintensity ///
c.smokeyrs##c.smokeyrs ib(last).exercise ib(last).active c.wt71##c.wt71 

predict pd_qsmk, pr

* estimation of numerator of ip weights 
logit qsmk 
predict pn_qsmk, pr

gen sw_a=.
replace sw_a=pn_qsmk/pd_qsmk if qsmk==1
replace sw_a=(1-pn_qsmk)/(1-pd_qsmk) if qsmk==0

* Second, fit logistic model below
logistic death qsmk [pweight=sw_a], cluster(seqn) 

/*NOTE:Stata has two commands for logistic regression, logit and logistic. 
The former displays the coefficients while the latter displays the odds ratios. 
You can also obtain the odds ratios by using the logit command with the or option/




***************************************************************
*PROGRAM 12.6
*Assessing effect modification by sex using a marginal structural mean model
*Data from NHEFS
***************************************************************/
clear
use "P:\CIBook\NHEFSdata\nhefs.dta"

* provisionally ignore subjects with missing values 
foreach var of varlist qsmk sex race age school smokeintensity smokeyrs exercise active wt71 wt82 {
drop if `var'==.
}


tab sex

/* estimation of denominator of ip weights */
logit qsmk sex race c.age##c.age ib(last).education c.smokeintensity##c.smokeintensity ///
c.smokeyrs##c.smokeyrs ib(last).exercise ib(last).active c.wt71##c.wt71 

predict pd_qsmk, pr

/* estimation of numerator of ip weights */
logit qsmk sex
predict pn_qsmk, pr

gen sw_a=.
replace sw_a=pn_qsmk/pd_qsmk if qsmk==1
replace sw_a=(1-pn_qsmk)/(1-pd_qsmk) if qsmk==0

summarize sw_a
regress wt82_71 qsmk##sex [pw=sw_a], cluster(seqn)



/***************************************************************
PROGRAM 12.7
Estimating IP weights to adjust for selection bias due to censoring
Data from NHEFS
***************************************************************/

clear
use "P:\CIBook\NHEFSdata\nhefs.dta"

*Analysis restricted to N=1629. Missing values of wt82 are retained while the others are dropped 

foreach var of varlist qsmk sex  race school active exercise age  wt71  smokeintensity smokeyrs {
drop if `var'==.
}


gen byte cens = (wt82 == .)

tab cens qsmk, column
bys cens: summarize wt71

/*Estimation of the denominator of ip weight for A*/
logit qsmk sex race c.age##c.age ib(last).education c.smokeintensity##c.smokeintensity ///
c.smokeyrs##c.smokeyrs ib(last).exercise ib(last).active c.wt71##c.wt71 

predict pd_qsmk, pr

/* estimation of numerator of ip weights for A */
logit qsmk
predict pn_qsmk, pr


/* estimation of denominator of ip weights for C */
logit cens qsmk sex race c.age##c.age ib(last).education c.smokeintensity##c.smokeintensity ///
c.smokeyrs##c.smokeyrs ib(last).exercise ib(last).active c.wt71##c.wt71 

predict pd_cens, pr

/* estimation of numerator of ip weights for C */
logit cens qsmk

predict pn_cens, pr

/*Estimation of stabilized weight for A (sw_a)*/

gen sw_a=.
replace sw_a=pn_qsmk/pd_qsmk if qsmk==1
replace sw_a=(1-pn_qsmk)/(1-pd_qsmk) if qsmk==0

/*Estimation of stabilized weight for C (sw_c)*/

/*NOTE: STATA, by default gives the probability of an event being 1. However, we are interested in Pr[C=0|A,L]
Observations with cens=1 only contribute to censoring models*/

gen sw_c=.
replace sw_c=(1-pn_cens)/(1-pd_cens) if cens==0


/*Estimation of Stabilized Censoring weight (sw)*/
gen sw=sw_a*sw_c
summarize sw
regress wt82_71 qsmk [pw=sw], cluster(seqn)
