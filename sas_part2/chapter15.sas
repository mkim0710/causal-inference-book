
/*******************************************************************
PROGRAM 15.1
Estimating the average causal effect within levels of confounders
under the assumption of effect-measure modification by smoking intensity ONLY
Data from NHEFS
********************************************************************/;
libname causinf "C:\dropbox\ci\data";
* ods html close; * use this option for output as text rather than html;

/* some preprocessing of the data */
data nhefs;
    set causinf.nhefs;
    cens= (wt82 eq .);
run;

/* regression on covariates, allowing for some effect modification  */
proc genmod data = nhefs;
    class exercise active education;
    model wt82_71 = qsmk sex race age age*age education
                smokeintensity smokeintensity*smokeintensity smokeyrs smokeyrs*smokeyrs
                exercise active wt71 wt71*wt71
                qsmk*smokeintensity;
    estimate 'smokeintensity = 5'  qsmk 1 qsmk*smokeintensity 5;
    estimate 'smokeintensity = 40' qsmk 1 qsmk*smokeintensity 40; 
run;
quit;

/* regression on covariates, not allowing for effect modification  */
proc genmod data = nhefs;
    class exercise active education;
    model wt82_71 = qsmk sex race age age*age education
                smokeintensity smokeintensity*smokeintensity smokeyrs smokeyrs*smokeyrs
                exercise active wt71 wt71*wt71;       
run;
quit;


/*******************************************************************
PROGRAM 15.2
Estimating and plotting the propensity score 
Data from NHEFS
********************************************************************/;
proc logistic data= nhefs descending;
    ods exclude ClassLevelInfo Type3 Association FitStatistics GlobalTests;
    class exercise active education;
    model qsmk = sex race age age*age education
                smokeintensity smokeintensity*smokeintensity smokeyrs smokeyrs*smokeyrs
                exercise active wt71 wt71*wt71;
    output out= est_prob p= ps;
run;

proc univariate data= est_prob;
	class qsmk;
	id seqn qsmk;
    var ps;
run;


/* plotting the estimated propensity score */
data tmp1 tmp2 ;
	set est_prob (keep = qsmk ps) ;
	ps = round(ps,0.05);
	if qsmk = 0 then output tmp1 ;
	if qsmk = 1 then output tmp2 ;
run;

proc freq data = tmp1 noprint  ;
	table ps / out = qsmk0 (keep = ps count ) ;
run;

proc freq data = tmp2 noprint ;
	table ps / out = qsmk1 (keep = ps count );
run;

data tmp3 ; 
	merge qsmk0 (in = a rename = (count = count0) ) qsmk1  (in = b rename = (count = count1));
	by ps ;
	label ps="Estimated Propensity Score"
    	  count0="No. Subjects"
    	  count1="No. Subjects" 
      	;
	if count0 = . then count0 = 0 ;
	if count1 = . then count1 = 0 ;
	count1 = -1 * count1 ;
	count2 = -1 *count1 ;
run;
 
proc format;
   picture positive 
     low-<0='000,000'
     0<-high='000,000';
run; 
ods graphics / reset border=off width=600px height=400px imagefmt=pdf  ;
title ;	footnote ;

data sganno;
	function ="text" ; label ="A=0" ; textcolor ="black" ;  x1 = 95 ; y1 = 97 ; size=15 ; 
	x1space="wallpercent" ; y1space = "wallpercent" ; output ;
	function ="text" ; label ="A=1" ; textcolor ="black" ;  x1 = 95 ; y1 = 3  ; size=15 ; 
	x1space="wallpercent" ; y1space = "wallpercent" ; output ;
run ;
 
proc sgplot data=tmp3 sganno=sganno noautolegend ;
    format count0 count1  positive.  ;
    vbar ps / response=count0   legendlabel=" "   datalabel datalabelattrs=(size=10)   name='a0' nofill ;
    vbar ps / response=count1   legendlabel=" "   datalabel datalabelattrs=(size=10)   name='a1'   fill fillattrs=(color=gray) ;
    xaxis labelattrs=(size=15) valueattrs=(size = 10) min = 0 max = 1.0;
    yaxis labelattrs=(size=15) valueattrs=(size = 10) ;
run;
 

/************************************************************************
PROGRAM 15.3
Stratification and outcome regression using deciles of the propensity score 
Data from NHEFS
*************************************************************************/;

/* calculation of deciles */
proc rank data=est_prob out=est_prob groups=10;
	var ps;
	ranks ps_dec;
proc sort; by ps_dec; run;

proc means data=est_prob;
	class ps_dec qsmk;
	var ps;
run;

/* stratification on PS deciles, allowing for effect modification */
proc ttest data= est_prob plots= none;
	by ps_dec;
	class qsmk;
	var wt82_71;
run;

/* regression on PS deciles, not allowing for effect modification */
proc glm data= est_prob;
	class ps_dec;
	model wt82_71= qsmk ps_dec /clparm solution;
run;


/*******************************************************************
PROGRAM 15.4
Standardization and outcome regression using the propensity score
Data from NHEFS
********************************************************************/;

/* regression on propensity score, not allowing for effect modification  */
proc glm data= est_prob;
	model wt82_71= qsmk ps/clparm;
run;
quit;


/* standardization by propensity score, agnostic regarding effect modification */;

%let nboot = 500 ; * number of bootstrap samples to run - to test the program use 50 or so ;

options nonotes nocenter ;

proc format ;
   value interv -1= "Observed"
                 0= "No treatment"
                 1= "Treatment"
                 2= "Treatment - No treatment"               
       ; 
run;

/* create a dataset with 3 copies of each subject */  
data onesample ;
  set est_prob end = _end_  ;
  label interv= "Intervention"; 
  retain _id ;
  if _n_ = 1 then _id = 0;
  _id = _id + 1 ;
  if _end_ then do ;
     call symput("nids",trim(left(_id)));
  end;
   
  interv = -1 ;    /* 1st copy: equal to original one */
    output ; 
  interv = 0 ;     /* 2nd copy: treatment set to 0, outcome to missing */
    qsmk = 0 ;
    wt82_71 = . ;
    output ;  
  interv = 1 ;     /* 3rd copy: treatment set to 1, outcome to missing*/
    qsmk = 1 ;
    wt82_71 = . ;
    output ;    
run;

* creation of bootstrap samples;
data ids ;
   do bsample = 1 to &nboot;
       do _id = 1 to &nids ;
           output ;
       end;
   end;
run;

proc surveyselect data= ids 
    method = urs
    n= &nids
    seed = 1232  
    out = _idsamples (keep = bsample _id  numberhits  ) 
    outall  noprint  ;
    strata bsample ;
run;

data _idsamples ; 
	set _idsamples ;
    if bsample = 0 then numberhits = 1 ;
run;
 
 proc sql ;
   	create table allbsamples as
   	select * from _idsamples full join onesample 
   	on _idsamples._id = onesample._id 
   	order by bsample ,_id 
  	;
quit;

* linear model to estimate mean outcome conditional on treatment and confounders;
* parameters are estimated using original observations only (interv= -1) ;
* parameter estimates are used to predict mean outcome for observations with set treatment (interv=0 and interv=1);
proc genmod data =  onesample ;
   	class exercise active education;
	model wt82_71 = qsmk ps;
	output out = predicted_mean0 p = meanY ; 
run;

data predicted_mean0 ;
	set predicted_mean0 ;
	bsample = 0 ;
	numberhits = 1 ;
run;

ods listing select none ;
proc genmod data = allbsamples;
    class exercise active education;
	model wt82_71 = qsmk ps;
	output out = predicted_mean1 p = meanY ;
    freq numberhits ;
    by bsample ;
run;
ods listing ;

data predicted_mean ;
	set predicted_mean0 predicted_mean1;
	by bsample ;
run;

proc sort data = predicted_mean ;
	by bsample interv ;
run;

* estimate mean outcome in each of the groups interv=-1, interv=0, and interv=1;
proc means data = predicted_mean mean noprint ; 
  	var meanY ;  
  	by bsample interv ;
  	freq numberhits ;
  	output out = results (keep = bsample interv mean ) mean = mean ;
run;

proc transpose data = results out = for_diff (keep = bsample col2 col3);
	var mean ;
	by bsample ;
run;

data for_diff ;
	set for_diff ;
	mean = col3 - col2 ;
	interv = 2 ;
	keep bsample interv mean ;
run;

proc means data = for_diff (where = (bsample > 0)) noprint ;
	var mean ;
	by interv ;
	output out=diffstd (keep = interv std) std = std ;
run;

* for bootstrap find std of each mean for bsample > 0 ;
proc means data = results (where = (bsample > 0)) noprint ;
	class interv ;
	var mean;
	types interv ;
	output out = stderrs (keep = interv std ) std = std ;
run;

data sample0;
	set results (where = (bsample = 0)) for_diff (where = (bsample = 0)) ;
	drop bsample ;
run;

data stderrs ;
	set stderrs diffstd ;
run;

data finalres;
	merge sample0 stderrs ;
	by interv ;
	lb = mean - 1.96 * std ;
	ub = mean + 1.96 * std ;
	label lb="95% Lower bound"
    	  ub="95% Upper bound"
      	  std="Standard Error"
     ;
run;

proc print data = finalres label noobs ;
	title1 "Parametric g-formula";
	title2 "Bootstrap results using &nboot samples" ;
	format interv interv. ;
	var interv mean std lb ub ;
run;

