 /***************************************************************
PROGRAM 16.1
Estimating the average causal effect using the standard IV estimator
via the calculation of sample averages
Data from NHEFS
***************************************************************/

clear
use "P:\CIBook\NHEFSdata\nhefs.dta"

gen byte cens = (wt82 == .)

drop if education ==.

summarize price82

/* ignore subjects with missing outcome or missing instrument for simplicity*/
foreach var of varlist wt82 price82 qsmk {
drop if `var'==.
}
/**/

/*Create categorical instrument*/
gen byte highprice  = (price82 >= 1.5 & price82 < .)

/*Calculate P[Z|A=a]*/
tab highprice qsmk, row

/*Calculate P[Y|Z=z]*/
ttest wt82_71, by(highprice)


/***************************************************************
PROGRAM 16.2
Estimating the average causal effect using the standard IV estimator
via two-stage-least-squares regression
Data from NHEFS
***************************************************************/

/*ivregress fits the model in two stages: */
/*first model: qsmk = highprice*/
/*second model: wt82_71 = predicted_qsmk*/
ivregress 2sls wt82_71 (qsmk = highprice)




/***************************************************************************
PROGRAM 16.3
Estimating the average causal effect using the standard IV estimator 
via an additive marginal structural model
Data from NHEFS
***************************************************************************/

gen psi = 2.396
gen hspi = wt82_71 -psi*qsmk

logit highprice hspi


/***************************************************************************
PROGRAM 16.4
Estimating the average causal effect using the standard IV estimator
based on alternative proposed instruments
Data from NHEFS
***************************************************************************/

/*Instrument cut-point: 1.6*/
replace highprice = .
replace highprice = (price82 >=1.6 & price82 < .)

ivregress 2sls wt82_71 (qsmk = highprice)


/*Instrument cut-point: 1.7*/
replace highprice = .
replace highprice = (price82 >=1.7 & price82 < .)

ivregress 2sls wt82_71 (qsmk = highprice)


/*Instrument cut-point: 1.8*/
replace highprice = .
replace highprice = (price82 >=1.8 & price82 < .)

ivregress 2sls wt82_71 (qsmk = highprice)


/*Instrument cut-point: 1.9*/
replace highprice = .
replace highprice = (price82 >=1.9 & price82 < .)

ivregress 2sls wt82_71 (qsmk = highprice)



/***************************************************************************
PROGRAM 16.5
Estimating the average causal effect using the standard IV estimator
conditional on baseline covariates
Data from NHEFS
***************************************************************************/

replace highprice = .
replace highprice = (price82 >=1.5 & price82 < .)

ivregress 2sls wt82_71 sex race age smokeintensity smokeyrs exercise active wt7 (qsmk = highprice)
