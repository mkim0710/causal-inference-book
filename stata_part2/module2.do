/**Code: Module 2***/

/***Part 1: Data set-up***/
clear
use "C:\Users\EMURRAY\Dropbox\CIBook\NHEFSdata\nhefs.dta"

/*use "<path>\nhefs.dta"*/

gen byte cens = (wt82 == .)
gen byte older = (age > 50 & age < .)
label var older "50+ years"
recode school (0/8 = 1 "1. 8th grade or less") (9/11 = 2 "2. HS dropout")(12 = 3 "3. HS") ///
(13/15 = 4 "4. College dropout") (16/max = 5 "5. College or more") (missing = 6 "Unknown"), gen(education)

/* provisionally ignore subjects with missing values*/
*Analysis restricted to N=1566 with no missing values*/
foreach var of varlist school wt82 {
drop if `var'==.
}

/***************/


/***Part 2: Non-parametric outcome regression****/

/*Non-parametric estimation in women and men separately*/
by sex, sort: tabulate qsmk, summarize(wt82_71)

by sex, sort: regress wt82_71 qsmk

/*Non-parametric estimation in women and men with a saturated linear regression model*/
/*E[Y|A,L]=theta0 + theta1A + theta2L +theta3A*L*/
gen qsmksex=qsmk*sex
regress wt82_71 qsmk sex qsmksex
predict meanY
by qsmk sex, sort: summarize meanY
drop meanY

/***************/


/****Part 3: Parametric outcome regression***/

/*Parametric model: E[Y|A,L]=theta0 + theta1A + theta2L */
regress wt82_71 qsmk sex
predict meanY
by qsmk sex, sort: summarize meanY
drop meanY


/***************/


/***Part 4: Estimating the propensity score***/

/*Non-parametric estimation of the PS using sample proportions*/
gen older = (age>50 & age < . )
by sex, sort: tabulate older, summarize(qsmk)


/*Non-parametric estimation of the PS using a saturated logistic regression model*/
logit qsmk sex##older 
predict meanY
by sex older, sort: summarize meanY
drop meanY


/*Parametric estimation of the PS using a non-saturated logistic regression model*/
logit qsmk sex older 
predict meanY
by sex older, sort: summarize meanY
drop meanY



/***************/



/***Part 5: Propensity score analyses***/

/*Estimate the propenity score*/
logit qsmk sex age
predict ps, pr
by qsmk, sort: summarize ps

/*create categorical cPS in deciles from the continuous PS*/
xtile ps_dec = ps, nq(10)
by ps_dec qsmk, sort: summarize ps


/*Stratify subjects by cPS category, and*/
/*Estimate the effect in each category*/
by ps_dec: ttest wt82_71, by (qsmk)

/*Regression on PS categories*/
regress wt82_71 qsmk i.ps_dec

/*Regression on continuous PS*/
regress wt82_71 qsmk ps












