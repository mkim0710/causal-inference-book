
/*******************************************************************
PROGRAM 16.1
Estimating the average causal effect using the standard IV estimator 
via the calculation of sample averages 
Data from NHEFS
********************************************************************/;
libname causinf "C:\dropbox\ci\data";
* ods html close; * use this option for output as text rather than html;

/* some preprocessing of the data */
data nhefs;
    set causinf.nhefs;
    cens= (wt82 eq .);
run;
proc univariate; var price82; run;

data nhefs_iv;
	set nhefs;
	* for simplicity, ignore subjects with missing outcome or missing instrument;
	where wt82 ne . and price82 ne .;
	highprice =(price82 ge 1.5); 
run;

proc freq;
	table highprice * qsmk;
run;

proc ttest plots= none;
	class highprice;
	var wt82_71;
run;

/*******************************************************************
PROGRAM 16.2
Estimating the average causal effect using the standard IV estimator 
via two-stage-least-squares regression
Data from NHEFS
********************************************************************/;

proc syslin data= nhefs_iv 2sls first;
   endogenous  qsmk; /* treatment variable */
   instruments highprice;
   model wt82_71= qsmk;  /* qsmk is actually replaced by its predicted value*/
run;


/*****************************************************************
PROGRAM 16.3
Estimating the average causal effect using the standard IV estimator 
via an additive marginal structural model
Data from NHEFS
******************************************************************/;


/*********************************************************/
/* G-estimation: Checking one possible value of psi      */
/* See Chapter 14 for program that checks several values */
/* and computes 95% confidence intervals                 */
/*********************************************************/
data nhefs_ivgest;
    set nhefs_iv;
    psi= 2.396;
    Hpsi= wt82_71-psi*qsmk;
run;

proc genmod data= nhefs_ivgest descending;
    ods exclude ClassLevels ParmInfo ModelInfo;  
    class seqn exercise active education;
    model highprice = 
				/* covariates can be added here */ 
				/* sex race age age*age education
                smokeintensity smokeintensity*smokeintensity smokeyrs smokeyrs*smokeyrs
                exercise active wt71 wt71*wt71 */ 
                Hpsi/ dist= bin link= logit;
    repeated subject=seqn / type=ind;
run;
quit;


/*******************************************************************
PROGRAM 16.4
Estimating the average causal effect using the standard IV estimator 
based on alternative proposed instruments
Data from NHEFS
********************************************************************/;

data nhefs_iv;
	set nhefs;
	where wt82 ne . and price82 ne .;
	highprice =(price82 ge 1.6); 
run;

proc syslin data= nhefs_iv 2sls first;
   endogenous  qsmk;
   instruments highprice;
   model wt82_71= qsmk; 
run;

data nhefs_iv;
	set nhefs;
	where wt82 ne . and price82 ne .;
	highprice =(price82 ge 1.7); 
run;

proc syslin data= nhefs_iv 2sls first;
   endogenous  qsmk;
   instruments highprice;
   model wt82_71= qsmk; 
run;

data nhefs_iv;
	set nhefs;
	where wt82 ne . and price82 ne .;
	highprice =(price82 ge 1.8); 
run;

proc syslin data= nhefs_iv 2sls first;
   endogenous  qsmk;
   instruments highprice;
   model wt82_71= qsmk; 
run;


data nhefs_iv;
	set nhefs;
	where wt82 ne . and price82 ne .;
	highprice =(price82 ge 1.9); 
run;

proc syslin data= nhefs_iv 2sls first;
   endogenous  qsmk;
   instruments highprice;
   model wt82_71= qsmk; 
run;



/*******************************************************************
PROGRAM 16.5
Estimating the average causal effect using the standard IV estimator 
conditional on baseline covariates
Data from NHEFS
********************************************************************/;

data nhefs_iv;
	set nhefs;
	where wt82 ne . and price82 ne .;
	highprice =(price82 ge 1.5); 
run;


proc syslin data= nhefs_iv 2sls first;
   endogenous  qsmk;
   instruments highprice sex race age smokeintensity smokeyrs
                exercise active wt71; /* proposed instrument plus baseline covariates */
   model wt82_71= qsmk sex race age smokeintensity smokeyrs
                exercise active wt71; 
run;

