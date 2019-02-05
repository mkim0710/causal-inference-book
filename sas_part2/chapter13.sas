
/***************************************************************
PROGRAM 13.1
Estimating the mean outcome within levels of treatment and confounders
Data from NHEFS
***************************************************************/;
libname causinf "C:\dropbox\ci\data";

/* some preprocessing of the data */
data nhefs;
	set causinf.nhefs;
	cens= (wt82 eq .);
run;

proc genmod data = nhefs;
	class exercise active education;
	model wt82_71 = qsmk sex race age age*age education
				smokeintensity smokeintensity*smokeintensity smokeyrs smokeyrs*smokeyrs
				exercise active wt71 wt71*wt71
				qsmk*smokeintensity;
    output out = predicted_mean p = meanY ;
run;
quit;

proc print; 
	where seqn= 24770;
	id seqn;
	var meanY qsmk sex race age education smokeintensity smokeyrs exercise active wt71;
run;

proc univariate;
	where cens=0; /* for comparison with observed outcome */
	var meanY wt82_71;
run;



/***************************************************************
PROGRAM 13.2
Standardizing the mean outcome to the baseline confounders
Data from Table 2.2
***************************************************************/;
data table22;
	input id$ L A Y interv;
	cards;
Rheia      0 0 0 -1
Kronos     0 0 1 -1
Demeter    0 0 0 -1
Hades      0 0 0 -1
Hestia     0 1 0 -1
Poseidon   0 1 0 -1
Hera       0 1 0 -1
Zeus       0 1 1 -1
Artemis    1 0 1 -1
Apollo     1 0 1 -1
Leto       1 0 0 -1
Ares       1 1 1 -1
Athena     1 1 1 -1
Hephaestus 1 1 1 -1
Aphrodite  1 1 1 -1
Cyclope    1 1 1 -1
Persephone 1 1 1 -1
Hermes     1 1 0 -1
Hebe       1 1 0 -1
Dionysus   1 1 0 -1
Rheia      0 0 . 0
Kronos     0 0 . 0
Demeter    0 0 . 0
Hades      0 0 . 0
Hestia     0 0 . 0
Poseidon   0 0 . 0
Hera       0 0 . 0
Zeus       0 0 . 0
Artemis    1 0 . 0
Apollo     1 0 . 0
Leto       1 0 . 0
Ares       1 0 . 0
Athena     1 0 . 0
Hephaestus 1 0 . 0
Aphrodite  1 0 . 0
Cyclope    1 0 . 0
Persephone 1 0 . 0
Hermes     1 0 . 0
Hebe       1 0 . 0
Dionysus   1 0 . 0
Rheia      0 1 . 1
Kronos     0 1 . 1
Demeter    0 1 . 1
Hades      0 1 . 1
Hestia     0 1 . 1
Poseidon   0 1 . 1
Hera       0 1 . 1
Zeus       0 1 . 1
Artemis    1 1 . 1
Apollo     1 1 . 1
Leto       1 1 . 1
Ares       1 1 . 1
Athena     1 1 . 1
Hephaestus 1 1 . 1
Aphrodite  1 1 . 1
Cyclope    1 1 . 1
Persephone 1 1 . 1
Hermes     1 1 . 1
Hebe       1 1 . 1
Dionysus   1 1 . 1
;
run;

proc genmod data = table22;
	class L;
	model Y = A L A*L;
    output out = predicted_mean p = meanY ;
run;

proc means data = predicted_mean mean noprint;
  class interv ;
  var meanY ;
  types interv ;
  output out = results (keep = interv mean ) mean = mean ;
run;

proc print data = results noobs label ;
  title "Standardized means";
  var interv mean ;
run;



/***************************************************************
PROGRAM 13.3
Standardizing the mean outcome to the baseline confounders
Data from NHEFS
***************************************************************/;


/* create a dataset with 3 copies of each subject */  
data onesample ;
  set nhefs ;
  label interv= "Intervention"; 
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


* linear model to estimate mean outcome conditional on treatment and confounders;
* parameters are estimated using original observations only (interv= -1) ;
* parameter estimates are used to predict mean outcome for observations with 
  treatment set to 0 (interv=0) and to 1 (innterv=1);
proc genmod data = onesample;
	class exercise active education;
	model wt82_71 = qsmk sex race age age*age education
				smokeintensity smokeintensity*smokeintensity smokeyrs smokeyrs*smokeyrs
				exercise active wt71 wt71*wt71;
    output out = predicted_mean p = meanY ;
run;

* estimate mean outcome in each of the groups interv=0, and interv=1;
* this mean outcome is a weighted average of the mean outcomes in each combination 
	of values of treatment and confounders, that is, the standardized outcome;
proc means data = predicted_mean mean noprint;
  class interv ;
  var meanY ;
  types interv ;
  output out = results (keep = interv mean ) mean = mean ;
run;

proc print data = results noobs label ;
  title "Parametric g-formula";
  var interv mean ;
run;



/*******************************************************************************
PROGRAM 13.4
Computing the 95% confidence interval of the standardized means and their difference
Data from NHEFS
********************************************************************************/;


%let nboot = 1000 ; * number of bootstrap samples to run - to test the program use 50 or so ;

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
  set nhefs end = _end_  ;
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
	model wt82_71 = qsmk sex race age age*age education
				smokeintensity smokeintensity*smokeintensity smokeyrs smokeyrs*smokeyrs
				exercise active wt71 wt71*wt71;
	output out = predicted_mean0 p = meanY ; 
run;

data predicted_mean0 ;
	set predicted_mean0 ;
	bsample = 0 ;
	numberhits = 1 ;
run;

ods listing select none ;
proc genmod data =  allbsamples ;
    class exercise active education;
	model wt82_71 = qsmk sex race age age*age education
				smokeintensity smokeintensity*smokeintensity smokeyrs smokeyrs*smokeyrs
				exercise active wt71 wt71*wt71;
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
