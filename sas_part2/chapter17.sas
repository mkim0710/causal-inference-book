/***************************************************************
PROGRAM 17.1
Nonparametric estimation of survival curves
Data from NHEFS
***************************************************************/;
libname causinf "C:\dropbox\ci\data";
ods graphics on ;

/* some preprocessing of the data */
data nhefs;
    set causinf.nhefs;
    if death=0 then survtime=120; 
	else if death=1 then survtime= (yrdth-83)*12 + modth; * yrdth ranges from 83 to 92;
run;

proc freq; table death*qsmk; run;
proc univariate; where death=1; var survtime; run;

proc lifetest data=nhefs plots= (s (atrisk=0 to 120 by 12 atrisktick),h) intervals= (0 to 120 by 1);
   time SurvTime*death(0);
   strata qsmk /test=logrank ;
   label survtime='Months of follow-up';
run;

/***************************************************************
PROGRAM 17.2
Parametric estimation of survival curves via hazards model
Data from NHEFS
***************************************************************/;

/* creation of person-month data */
data nhefs_surv;
	length seqn 8. time 8. event 8.;
   	set nhefs;
	do time= 0 to (survtime-1);
		event= (death=1 and time=survtime-1);
		timesq= time*time;
		output;
	end;
run;

/* fit of hazards model */
proc logistic data= nhefs_surv outmodel=hazards_model;
	ods exclude ClassLevelInfo Type3 Association FitStatistics GlobalTests;
	model event = qsmk qsmk*time qsmk*timesq time timesq;
	output out=est_prob0 p=p_noevent;
run;

/* creation of dataset with all time points under each treatment level */
data for_hazards ;
	do time = 0 to 119;
	    timesq = time * time ;
		qsmk = 0 ;
		output ;
		qsmk = 1 ;
		output ;
	end;
run;

/* assignment of estimated (1-hazard) to each person-month */
proc logistic inmodel=hazards_model;
	score data=for_hazards out=pred  (keep = time qsmk P_0 rename = (P_0=p_noevent)) ;
run;
proc sort data = pred ;	by qsmk time ; run;

/* computation of survival for each person-month */
data pred ;
	set pred ;
	by qsmk ;
	retain surv ;
	if first.qsmk then surv = 1 ;
	surv = surv * p_noevent ;
run;

/* some data management to plot estimated survival curves */
data qsmk0 qsmk1 ;
	set pred ;
	if qsmk = 0 then output qsmk0 ;
	if qsmk = 1 then output qsmk1 ;
	keep time surv;
run;

data hazards_graph ;
	merge qsmk0 (rename = (surv = surv0)) qsmk1 (rename = (surv = surv1));
	survdiff= surv1 - surv0;
run;

proc sgplot data = hazards_graph ;
   series x = time y = surv0 /legendlabel='A = 0';
   series x = time y = surv1 /legendlabel='A = 1';
   title 'Survival from hazards model';
   yaxis label = 'Survival '  values = (0.6 to 1 by 0.2) ;
   xaxis label = 'Months '  values = (0 to 120 by 12) ;
   keylegend/noborder title='A: ';
run;
  

/***************************************************************
PROGRAM 17.3
Estimation of survival curves via IP weighted hazards model
Data from NHEFS
***************************************************************/;

/* estimation of denominator of IP weights */
proc logistic data=nhefs descending;
	ods exclude ClassLevelInfo Type3 Association FitStatistics GlobalTests;
	class exercise active education;
	model qsmk = sex race age age*age education
				smokeintensity smokeintensity*smokeintensity smokeyrs smokeyrs*smokeyrs
				exercise active wt71 wt71*wt71;
	output out=est_prob_d p=pd_qsmk;
run;
proc sort; by seqn; run;

/* estimation of numerator of ip weights */
proc logistic data=nhefs descending;
	ods exclude ClassLevelInfo Type3 Association FitStatistics GlobalTests Oddsratios;
	model qsmk = ;
	output out=est_prob_n (keep= seqn pn_qsmk) p=pn_qsmk;
run;
proc sort; by seqn; run;

/* computation of estimated weights */
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

/* creation of person-month data */
data nhefs_ipw;
	length seqn 8. time 8. event 8.;
   	set nhefs_sw;
	do time= 0 to (survtime-1);
		event= (death=1 and time=survtime-1);
		timesq= time*time;
		output;
	end;
run;

/* fit of weighted hazards model */
proc logistic data= nhefs_ipw outmodel=ipw_model;
	ods exclude ClassLevelInfo Type3 Association FitStatistics GlobalTests;
	weight sw_a;
	model event = qsmk qsmk*time qsmk*timesq time timesq;
	output out=est_prob0 p=p_noevent;
run;

/* creation of survival curves */
data for_ipw ;
	do time = 0 to 119;
	    timesq = time * time ;
		qsmk = 0 ;
		output ;
		qsmk = 1 ;
		output ;
	end;
run;

proc logistic inmodel=ipw_model ;
	score  data=for_ipw out= ipw_pred (keep = time qsmk P_0 rename = (P_0=p_noevent)) ;
run;
proc sort data = ipw_pred ;	by qsmk time ; run;

data ipw_pred ;
	set ipw_pred ;
	by qsmk ;
	retain surv ;
	if first.qsmk then surv = 1 ;
	surv = surv * p_noevent ;
run;

data qsmk0 qsmk1 ;
	set ipw_pred ;
	if qsmk = 0 then output qsmk0 ;
	if qsmk = 1 then output qsmk1 ;
	keep time surv;
run;

data ip_graph ;
	merge qsmk0 (rename = (surv = surv0)) qsmk1 (rename = (surv = surv1));
	survdiff= surv1 - surv0;
run;
proc print data= ip_graph; id time; run;

proc sgplot data = ip_graph ;
   series x = time y = surv0 /legendlabel='A = 0';
   series x = time y = surv1 /legendlabel='A = 1';
   title 'Survival from IP weighted hazards model';
   yaxis label = 'Survival '  values = (0.6 to 1 by 0.2) ;
   xaxis label = 'Months '  values = (0 to 120 by 12) ;
   keylegend/noborder title=' ';
run;
title;

/***************************************************************
PROGRAM 17.4
Estimating of survival curves via g-formula
Data from NHEFS
***************************************************************/;

/* fit of hazards model with covariates */
proc logistic data= nhefs_surv outmodel = gf_model;
	ods exclude ClassLevelInfo Type3 Association FitStatistics GlobalTests;
	class exercise active education;
	model event = qsmk qsmk*time qsmk*timesq time timesq
				sex race age age*age education
				smokeintensity smokeintensity*smokeintensity smkintensity82_71
				smokeyrs smokeyrs*smokeyrs
				exercise active wt71 wt71*wt71;
	output out=est_prob p=p_noevent;
run;

/* creation of dataset with all time points for each individual under each treatment level */
data for_gf ;
	set nhefs_surv ( where = (time = 0));
	do time = 0 to 119;
	    timesq = time*time ;
		qsmk = 0 ;
		output ;
		qsmk = 1 ;
		output ;
	end;
run;
 
proc logistic inmodel=gf_model ;
	score  data=for_gf out=gf_pred  (keep = seqn time qsmk P_0 rename = (P_0=p_noevent)) ;
run;
proc sort data = gf_pred ; by seqn qsmk time ; run;

data gf_pred ;
	set gf_pred ;
	by seqn qsmk ;
	retain surv ;
	if first.qsmk then surv = 1 ;
	surv = surv * p_noevent ;
run;

proc means data = gf_pred noprint  ;
	class qsmk time ;
	var surv ;
	types qsmk*time ;
	output out= mysurv_gf (keep  = qsmk time surv) mean=surv ;
run;

data qsmk0_gf qsmk1_gf ;
	set mysurv_gf ;
	if qsmk = 0 then output qsmk0_gf ;
	if qsmk = 1 then output qsmk1_gf ;
	keep time surv;
run;

data gf_graph ;
	merge qsmk0_gf (rename = (surv = surv0)) qsmk1_gf (rename = (surv = surv1));
	survdiff= surv1 - surv0;
run;
proc print data= gf_graph; id time; run;

proc sgplot data = gf_graph ;
   series x = time y = surv0 /legendlabel='A = 0';
   series x = time y = surv1 /legendlabel='A = 1';
   title 'Survival from g-formula';
   yaxis label = 'Survival '  values = (0.6 to 1 by 0.2) ;
   xaxis label = 'Months '  values = (0 to 120 by 12) ;
   keylegend/noborder title='A: ';
run;
title;

/***************************************************************
PROGRAM 17.5
Estimating of median survival time ratio via a structural nested AFT model
Data from NHEFS
***************************************************************/;

/* some preprocessing of the data */
data nhefs;
    set causinf.nhefs_book;
    if education ne .;
	if death=0 then survtime= .; 
	else if death=1 then survtime= (yrdth-83)*12 + modth; * yrdth ranges from 83 to 92;
run;

/* model to estimate E[A|L] */
proc logistic data=nhefs descending;
	ods exclude ClassLevelInfo Type3 Association FitStatistics GlobalTests;
	class exercise active education;
	model qsmk = sex race age age*age education
				smokeintensity smokeintensity*smokeintensity smokeyrs smokeyrs*smokeyrs
				exercise active wt71 wt71*wt71;
	output out=_modelA_ (keep = seqn  p_qsmk) p=p_qsmk;
run;
 
data  _uno_  ;
    merge nhefs _modelA_  end= _end_;
    by seqn;
    if  . < survtime < 120  then _events_used + 1 ; 
    if _end_ then call symput('events_used',left(_events_used));                                     
    if .< survtime <= 120 ; * select only those with observed death time;
run;

proc iml; 

    /* Compute the estimating function function that needs to be minimized */
    start sumeef(Psi);
    	load qsmk p_qsmk Amax Amin survtime;
        rows = nrow(survtime);
        Smat = J(rows,1,0);
               
		/* creation of delta indicator */
        do i = 1 to rows ;
		    T = survtime[i];
			delta = 0 ;
			if psi >= 0 then do;
			      if qsmk[i] = 0 then delta = 1 ;
				  if qsmk[i] = 1 & psi <= log(120/T) then delta = 1 ;
			end;
            else if psi < 0 then do ;
                  if qsmk[i] = 1 then delta = 1 ;
                  if qsmk[i] = 0 & psi > log(T/120) then delta = 1 ;
            end;             
            EA = p_qsmk[i];             
            Smat[i] =  delta * (qsmk[i] - EA) ;                                                                         
        end;

        Sval = Smat[+,];
        Save = Sval/rows;
        Smat = Smat - repeat(Save,rows);

        /* covariance */
		Sigma = 0;
        sigma = t(smat)*smat ;
	    if sigma = 0 then sigma = 1.0e-16;
	    estimeq = (Sval)*Inv(Sigma)*Sval`;
                           
        return(estimeq);

	finish sumeef;

    use _uno_;
    read all var{ qsmk  p_qsmk   survtime};    
    close _uno_;
     
    Amin = min(qsmk);
    Amax = max(qsmk);
    store Amin Amax  ;                           
    store qsmk  p_qsmk  survtime;                                                                                                   

	con = { -1  , 1 } ; /* optimization bounds */         
	optn = {0 1};       /* minimization print (all: 5, nothing: 0 */            
	Psi0 = 0;           /* starting value of psi for minimization procedure */
	call nlpnms(rc,xr,"SUMEEF",psi0,optn,con); * Nelder-Mead Simplex optimization;
	psi1= xr[1]; 
	objfunc = sumeef(xr);

	create _kk_ var {objfunc psi1 };
	append var {objfunc psi1 };
	close _kk_;
   
    start for_conf(x);
        tmp = sumeef(x) -3.84;
        return (tmp);
    finish for_conf;

 	* use simple bisection method to find estimates of lower and upper 95% confidence bounds ;
    increm = 0.1;
    psilow = xr;
    testlow = objfunc;  
    countlow = 0;
    if objfunc < 3.84 then do ;
		* find estimate of where sumeef(x) > 3.84 ;
	    do while ( testlow < 3.84 & countlow < 100 );
	        psilow = psilow - increm;
	        testlow = sumeef(psilow);
	        countlow = countlow + 1; 
	    end;
	   
	 	/* upper bound of 95% CI */
	    psihigh = xr;
	    testhigh = objfunc;
	    counthigh = 0;
	    do while (testhigh < 3.84 & counthigh < 100 );
	        psihigh = psihigh + increm;
	        testhigh = sumeef(psihigh);
	        counthigh = counthigh + 1;
	    end;
	
		/* better estimate using bisection method */
		if (testhigh > 3.84) & (testlow > 3.84) then do;
           
			/* bisection method */
			left = xr;
			fleft = objfunc - 3.84;
			right = psihigh;
			fright = testhigh - 3.84 ;
			middle = (left + right)/2.0;
			fmiddle = for_conf(middle);
			count = 0;
			diff = right - left;

			do until(abs(fmiddle) < 0.0001 | diff < 0.0001 | count > 100 );
				test = fmiddle * fleft ;
				if test < 0 then do;
					right = middle;
					fright = fmiddle;
				end;
				else do;
					left = middle;
					fleft = fmiddle;
				end;
				middle = (left + right)/2.0;
				fmiddle = for_conf(middle);
				count = count + 1;        
				diff = right - left;
			end;

			psi_high = middle;
			objfunc_high = fmiddle + 3.84;
		 
			/* lower bound of 95% CI*/
			left = psilow;
			fleft = testlow - 3.84;
			right = xr;
			fright = objfunc - 3.84 ;
			middle = (left + right)/2.0;
			fmiddle = for_conf(middle);
			count = 0;
			diff = right - left ;
			do until ( abs(fmiddle) < 0.0001 | diff < 0.0001  | count > 100 );
				test = fmiddle * fleft ;
				if test < 0 then do;
					right = middle;
					fright = fmiddle;
				end;
				else do;
					left = middle;
					fleft = fmiddle;
				end;
				middle = (left + right)/2.0;
				fmiddle = for_conf(middle);
				diff = right - left ;
				count = count + 1;                            
			end;
			psi_low = middle;
			objfunc_low = fmiddle + 3.84; 
			psi = xr ;
			print psi psi_low psi_high;
            
       end;
	end;

quit;
title; footnote;

