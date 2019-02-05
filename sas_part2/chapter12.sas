/***************************************************************
PROGRAM 12.1
Descriptive statistics from NHEFS data
***************************************************************/;
libname causinf "C:\dropbox\ci\data";

/* some preprocessing of the data */
data nhefs;
	set causinf.nhefs;
	match=1;
	cens= (wt82 eq .);
run;

data nhefs_nmv;
	set nhefs;
	if wt82 ne .; * provisionally ignore subjects with missing values for weight in 1982;
run;

proc genmod data= nhefs_nmv; 
	model wt82_71= qsmk;
	estimate 'Smoking cessation' intercept 1 qsmk 1;
	estimate 'No smoking cessation' intercept 1 qsmk 0;
run;
quit;

/* Table */
proc means data= nhefs_nmv;
	class qsmk;
	var age wt71 smokeintensity smokeyrs;
run;
proc freq data= nhefs_nmv;
	table qsmk*(sex race education exercise active) / nopercent nocol;
run;


/***************************************************************
PROGRAM 12.2
Estimating IP weights
Data from NHEFS
***************************************************************/;

/* estimation of ip weights via a logistic model */
proc logistic data= nhefs_nmv descending;
	ods exclude ClassLevelInfo ModelAnova Association FitStatistics GlobalTests;
	class exercise active education;
	model qsmk = sex race age age*age education
				smokeintensity smokeintensity*smokeintensity smokeyrs smokeyrs*smokeyrs
				exercise active wt71 wt71*wt71;
	output out=est_prob p=p_qsmk;
run;

data nhefs_w;
	set est_prob;

	if qsmk=1 then w= 1/p_qsmk;
	else if qsmk=0 then w= 1/(1-p_qsmk);
	
run;

proc univariate data=nhefs_w;
	id seqn;
	var w;
run;

proc genmod data= nhefs_w; 
	class seqn;
	weight w;
	model wt82_71= qsmk;
	estimate 'Smoking cessation' intercept 1 qsmk 1;
	estimate 'No smoking cessation' intercept 1 qsmk 0;
	repeated subject=seqn / type=ind;
run;
quit;

/* no association between sex and qsmk in pseudo-population*/ 
proc freq data=nhefs_w;
	weight w;
	table sex*qsmk;
run;

/* "check" for positivity */
title 'White women';
proc freq data=nhefs_nmv;
	where sex=1 and race=0; 
	table age*qsmk;
run;
title;



/***************************************************************
PROGRAM 12.3
Estimating stabilized IP weights
Data from NHEFS
***************************************************************/;

/* estimation of denominator of ip weights */
proc logistic data=nhefs_nmv descending;
	ods exclude ClassLevelInfo ModelAnova Association FitStatistics GlobalTests;
	class exercise active education;
	model qsmk = sex race age age*age education
				smokeintensity smokeintensity*smokeintensity smokeyrs smokeyrs*smokeyrs
				exercise active wt71 wt71*wt71;
	output out=est_prob_d p=pd_qsmk;
run;
proc sort; by seqn; run;

/* estimation of numerator of ip weights */
proc logistic data=nhefs_nmv descending;
	ods exclude ClassLevelInfo ModelAnova Association FitStatistics GlobalTests Oddsratios;
	model qsmk = ;
	output out=est_prob_n (keep= seqn pn_qsmk) p=pn_qsmk;
run;
proc sort; by seqn; run;

data nhefs_sw;
	merge est_prob_d est_prob_n ;
	by seqn;

	if qsmk=1 then sw_a= pn_qsmk / pd_qsmk;
	else if qsmk=0 then sw_a= (1-pn_qsmk) / (1-pd_qsmk);

run;

proc univariate data=nhefs_sw;
	var sw_a;
	id seqn;
run;

proc genmod data= nhefs_sw; 
	class seqn;
	weight sw_a;
	model wt82_71= qsmk;
	estimate 'Smoking cessation' intercept 1 qsmk 1;
	estimate 'No smoking cessation' intercept 1 qsmk 0;
	repeated subject=seqn / type=ind;
run;
quit;

/* no association between sex and qsmk in pseudo-population*/ 
proc freq data=nhefs_sw;
	weight sw_a;
	table sex*qsmk;
run;



/***************************************************************
PROGRAM 12.4
Estimating the parameters of a marginal structural mean model
with a continuous treatment
Data from NHEFS
***************************************************************/;

* Analysis restricted to subjects reporting <=25 cig/day at baseline;
data nhefs_nmv_s;
	set nhefs_nmv;
	if smokeintensity le 25;
run;

/* estimation of denominator of ip weights */
proc glm data= nhefs_nmv_s 
         outstat= ss_den(keep= _source_ _type_ df ss where=(_source_ in('ERROR') and _type_ in('ERROR')));;
	class exercise active education;
	model smkintensity82_71 = sex race age age*age education
				smokeintensity smokeintensity*smokeintensity smokeyrs smokeyrs*smokeyrs
				exercise active wt71 wt71*wt71 / solution;
	output out= temp_den p= pred;
run;

data sd_den;
	set ss_den;
	rootmse_n= sqrt(ss/df);
	match= 1;
	keep rootmse_n match;
proc print; run;

data est_dens_d;
        merge temp_den sd_den;
		by match;
        dens_den = pdf('NORMAL', smkintensity82_71, pred, rootmse_n);
proc sort; by seqn; run;

/* estimation of numerator of ip weights */
proc glm data= nhefs_nmv_s 
         outstat= ss_num(keep= _source_ _type_ df ss where=(_source_ in('ERROR') and _type_ in('ERROR')));;
	model smkintensity82_71 = / solution;
	output out= temp_num p= pred;
run;

data sd_num;
	set ss_num;
	rootmse_n= sqrt(ss/df);
	match= 1;
	keep rootmse_n match;
proc print; run;

data est_dens_n;
        merge temp_num sd_num;
		by match;
        dens_num = pdf('NORMAL', smkintensity82_71, pred, rootmse_n);
proc sort; by seqn; run;

data nhefs_sw_cont;
	merge est_dens_d est_dens_n ;
	by seqn;
	sw_a= dens_num / dens_den;
run;

proc univariate data=nhefs_sw_cont;
	var sw_a;
	id seqn;
run;

proc genmod data= nhefs_sw_cont; 
	class seqn;
	weight sw_a;
	model wt82_71= smkintensity82_71 smkintensity82_71*smkintensity82_71;
	estimate 'No change' intercept 1 smkintensity82_71 0;
	estimate 'Increase smoking by 20 cig/day' intercept 1 smkintensity82_71 20;
	repeated subject=seqn / type=ind;
run;
quit;



/***************************************************************
PROGRAM 12.5
Estimating the parameters of a marginal structural logistic model
Data from NHEFS
***************************************************************/;

proc freq data= nhefs_nmv;
	table qsmk*death;
run;

* First, estimation of stabilized weights sw_a (same as in PROGRAM 12.3);
* Second, fit logistic model below;
proc genmod data= nhefs_sw descending; 
	class seqn;
	weight sw_a;
	model death = qsmk / link=logit dist=bin;
	estimate 'Log OR' qsmk 1 /exp;
	repeated subject=seqn / type=ind;
run;
quit;



/***************************************************************
PROGRAM 12.6
Assessing effect modification by sex using a marginal structural mean model
Data from NHEFS
***************************************************************/;

proc freq data= nhefs_nmv;
	table sex;
run;

/* estimation of denominator of ip weights */
proc logistic data=nhefs_nmv descending;
	ods exclude ClassLevelInfo ModelAnova Association FitStatistics GlobalTests;
	class exercise active education;
	model qsmk = sex race age age*age education
				smokeintensity smokeintensity*smokeintensity smokeyrs smokeyrs*smokeyrs
				exercise active wt71 wt71*wt71;
	output out=est_prob_d p=pd_qsmk;
run;
proc sort; by seqn; run;

/* estimation of numerator of ip weights */
proc logistic data=nhefs_nmv descending;
	ods exclude ClassLevelInfo ModelAnova Association FitStatistics GlobalTests Oddsratios;
	model qsmk = sex;
	output out=est_prob_n (keep= seqn pn_qsmk) p=pn_qsmk;
run;
proc sort; by seqn; run;

data nhefs_sw;
	merge est_prob_d est_prob_n ;
	by seqn;

	if qsmk=1 then sw_a= pn_qsmk / pd_qsmk;
	else if qsmk=0 then sw_a= (1-pn_qsmk) / (1-pd_qsmk);

run;

proc univariate data=nhefs_sw;
	var sw_a;
	id seqn;
run;

proc genmod data= nhefs_sw; 
	class seqn sex;
	weight sw_a;
	model wt82_71= qsmk sex qsmk*sex;
	repeated subject=seqn / type=ind;
run;
quit;



/***************************************************************
PROGRAM 12.7
Estimating IP weights to adjust for selection bias due to censoring
Data from NHEFS
***************************************************************/;

proc freq data=nhefs;
	table qsmk*cens;
run;
proc means data= nhefs;
	class cens;
	var wt71;
run;

/* estimation of denominator of ip weights for A */
proc logistic data=nhefs descending;
	ods exclude ClassLevelInfo ModelAnova Association FitStatistics GlobalTests;
	class exercise active education;
	model qsmk = sex race age age*age education
				smokeintensity smokeintensity*smokeintensity smokeyrs smokeyrs*smokeyrs
				exercise active wt71 wt71*wt71;
	output out=est_prob_d_a p=pd_qsmk;
run;
proc sort; by seqn; run;

/* estimation of numerator of ip weights for A */
proc logistic data=nhefs descending;
	ods exclude ClassLevelInfo ModelAnova Association FitStatistics GlobalTests Oddsratios;
	model qsmk = ;
	output out=est_prob_n_a (keep= seqn pn_qsmk) p=pn_qsmk;
run;
proc sort; by seqn; run;


/* estimation of denominator of ip weights for C */
proc logistic data=nhefs;
	ods exclude ClassLevelInfo ModelAnova Association FitStatistics GlobalTests;
	class exercise active education;
	model cens = qsmk sex race age age*age education
				smokeintensity smokeintensity*smokeintensity smokeyrs smokeyrs*smokeyrs
				exercise active wt71 wt71*wt71;
	output out=est_prob_d_c (keep= seqn pd_cens) p=pd_cens;
run;
proc sort; by seqn; run;

/* estimation of numerator of ip weights for C */
proc logistic data=nhefs;
	ods exclude ClassLevelInfo ModelAnova Association FitStatistics GlobalTests Oddsratios;
	model cens = qsmk;
	output out=est_prob_n_c (keep= seqn pn_cens) p=pn_cens;
run;
proc sort; by seqn; run;


data nhefs_sw;
	merge est_prob_d_a est_prob_n_a est_prob_d_c est_prob_n_c;
	by seqn;
	if cens=0; * observations with cens=1 only contribute to censoring models;

	if qsmk=1 then sw_a= pn_qsmk / pd_qsmk;
	else if qsmk=0 then sw_a= (1-pn_qsmk) / (1-pd_qsmk);
	sw_c= pn_cens / pd_cens;
	sw= sw_a * sw_c;

run;

proc univariate data=nhefs_sw;
	var sw_a sw_c sw;
	id seqn;
run;

proc genmod data= nhefs_sw; 
	class seqn;
	weight sw;
	model wt82_71= qsmk;
	estimate 'Smoking cessation' intercept 1 qsmk 1;
	estimate 'No smoking cessation' intercept 1 qsmk 0;
	repeated subject=seqn / type=ind;
run;
quit;


