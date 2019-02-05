libname causinf   'C:\dropbox\ci\data';
ods graphics on ;
 
/***************************************************************
PROGRAM 17.3b
Estimation of survival curves via IP weighted hazards model
with bootstrapping for variance estimation
Data from NHEFS
***************************************************************/;

/* bootstrap part - stack bootstrap data sets and use same method as in Program 17.3 */
data nhefs_surv;
    set causinf.nhefs;
	length sample 3 ;
	if death=0 then survtime=120; 
	else if death=1 then survtime= (yrdth-83)*12 + modth; * yrdth ranges from 83 to 92;
    keep qsmk  seqn sex race age education smokeintensity smokeyrs exercise active wt71 survtime death smkintensity82_71;
run;
 

%macro bootstrap(nboot = 500);

	data justids;
		set nhefs_surv (keep = seqn);
	run;

	data justids_boot ;
		set justids ;
		do sample = 0 to &nboot ;
	 	   output ;
		end;
	run;

	proc sort data = justids_boot ;	by sample seqn ; run;

	proc surveyselect data= justids_boot 
 	     method = urs
  	     n= 1629
   	     seed = 1232  
    	 out = justids_boot (drop = selected samplingweight expectedhits)
     	 outall    
         noprint               ;
  		strata sample ;
  	run;

	data justids_boot ;
		set justids_boot ;
		if sample = 0 then numberhits = 1 ;
	run;

	proc sort data = justids_boot ;
		by sample seqn ;
	run;
 
	%do sample = 0 %to &nboot ;
	    /* estimation of denominator of ip weights */
		data nhefs_boot ;
			merge nhefs_surv justids_boot(where = (sample = &sample));
			by seqn ;
			if numberhits = 0 then delete ;
		run;

		proc logistic data=nhefs_boot descending noprint ;
			ods exclude ClassLevelInfo Type3 Association FitStatistics GlobalTests;
			class exercise active education;
			model qsmk = sex race age age*age education
						smokeintensity smokeintensity*smokeintensity smokeyrs smokeyrs*smokeyrs
						exercise active wt71 wt71*wt71;
			output out=est_prob_d p=pd_qsmk;
			freq numberhits ;
		run;
		proc sort; by sample seqn; run;

	/* estimation of numerator of ip weights */
		proc logistic data=nhefs_boot descending noprint;
			ods exclude ClassLevelInfo Type3 Association FitStatistics GlobalTests Oddsratios;
			model qsmk = ;
			output out=est_prob_n (keep=  sample seqn pn_qsmk) p=pn_qsmk;
			freq numberhits ; 
		run;

		proc sort; by sample seqn; run;

		data nhefs_sw_boot;
			merge est_prob_d est_prob_n ;
			by sample seqn;
			if qsmk=1 then sw_a= pn_qsmk / pd_qsmk;
			else if qsmk=0 then sw_a= (1-pn_qsmk) / (1-pd_qsmk);
		run;

	/* creation of person-month data */
		data nhefs_surv_boot;
			length seqn 8. time 3. event 3.;
		   	set nhefs_sw_boot;
			do time= 0 to (survtime -1);
				event= (death=1 and time=(survtime-1));
				timesq= time*time;
				output;
			end;
		run;

		proc logistic data= nhefs_surv_boot outmodel=ipw_model noprint ;
			ods exclude ClassLevelInfo Type3 Association FitStatistics GlobalTests;
			weight sw_a;
			model event = qsmk qsmk*time time timesq;
			*output out=est_prob0 p=p_noevent;
			freq numberhits ;
		run;

		%if &sample = 0 %then %do;
			data for_ipw_boot ;
		 		length sample qsmk time 3 ;
				numberhits = 1 ;
				do time = 0 to 119;
				    timesq = time * time ;
					qsmk = 0 ;
					output ;
					qsmk = 1 ;
					output ;
				end;
	 		run;
		%end;

		proc logistic inmodel=ipw_model  ;
			score  data=for_ipw_boot out=pred_boot  (keep =  numberhits time qsmk P_0 rename = (P_0=p_noevent)) ;
		run;

		proc sort data = pred_boot ; by  qsmk time ; run;

		data pred_boot ;
			set pred_boot ;
			by  qsmk ;
			retain surv ;
			if first.qsmk then surv = 1 ;
			surv = surv * p_noevent ;
		run;

		data qsmk0 qsmk1 ;
			set pred_boot ;
			if qsmk = 0 then output qsmk0 ;
			if qsmk = 1 then output qsmk1 ;
			keep  time surv;
		run;

		data ip_graphs ;
			merge qsmk0 ( rename = (surv = surv0)) qsmk1 ( rename = (surv = surv1));
			by  time ;
			survdiff = surv1 - surv0 ;
			sample = &sample ;
		run;

		%if &sample = 0 %then %do ;
			data ip_graphs_boot ;
				set ip_graphs ;
			run;

		%end;
		%else %do;
		    proc append base = ip_graphs_boot data=ip_graphs;
			run;
		%end;

		proc datasets library = work nolist ;
			delete ip_graphs qsmk0 qsmk1 pred_boot nhefs_boot nhefs_sw_boot est_prob_n est_prob_d ipw_model for_ip_boot ;
		quit;
 
/***************************************************************
PROGRAM 17.4b
Estimating of survival curves via g-formula
with bootstrapping for variance estimation
Data from NHEFS
***************************************************************/;
 
    proc logistic data= nhefs_surv_boot outmodel = gf_model noprint ;
		ods exclude ClassLevelInfo Type3 Association FitStatistics GlobalTests;
		class exercise active education;
		model event = qsmk qsmk*time qsmk*timesq time timesq
					sex race age age*age education
					smokeintensity smokeintensity*smokeintensity smkintensity82_71
					smokeyrs smokeyrs*smokeyrs
					exercise active wt71 wt71*wt71;
		freq numberhits ;
	run;

	data for_gf ;
		set nhefs_surv_boot(where = (time = 0));
		do time = 0 to 119 ;
		     timesq = time*time ;
			 qsmk = 0 ;
			 output ;
			 qsmk = 1 ;
			 output ;
		end;
	run;

	proc logistic inmodel=gf_model ;
		score  data=for_gf out=gf_pred  (keep =  numberhits seqn time qsmk P_0 rename = (P_0=p_noevent)) ;
	run;

	proc sort data = gf_pred ;
		by seqn qsmk time ;
	run;

	data gf_pred ;
		set gf_pred ;
		by  seqn qsmk ;
		retain surv ;
		if first.qsmk then surv = 1 ;
		surv = surv * p_noevent ;
	run;

	proc means data = gf_pred noprint  ;
		class qsmk time ;
		var surv ;
		types qsmk*time ;
		output out= mysurv_gf (keep  =  qsmk time surv) mean=surv ;
		freq numberhits ;
	run;

	data qsmk0_gf qsmk1_gf ;
		set mysurv_gf ;
		if qsmk = 0 then output qsmk0_gf ;
		if qsmk = 1 then output qsmk1_gf ;
		keep  time surv;
	run;

	data gf_graphs ;
		merge qsmk0_gf ( rename = (surv = surv0)) qsmk1_gf ( rename = (surv = surv1));
		by  time ;
		surv_diff = surv1 - surv0 ;
		sample = &sample ;
	run;


		%if &sample = 0 %then %do ;
			data gf_graphs_boot ;
				set gf_graphs ;
			run;

		%end;
		%else %do;
		    proc append base = gf_graphs_boot data=gf_graphs;
			run;
		%end;

		proc datasets library = work nolist ;
			delete gf_graphs qsmk0_gf qsmk1_gf pred_boot nhefs_surv_boot   gf_model for_gf gf_pred  ;
		quit;

	 %end ;


	data sample0 ;
		set ip_graphs_boot (where = (sample = 0));
	run;

	proc sort data = ip_graphs_boot ; by time ; run;

	proc univariate data= ip_graphs_boot(where = (sample > 0)) noprint;
		var survdiff ;
		by time ;
		output out=mysurv_pct pctlpre=diff_   pctlpts=2.5,97.5; 
	run; 

	data ip_graph ;
		merge sample0  mysurv_pct ;
		by time ;
	run;

	proc print data= ip_graph; id time; run;

	proc sgplot data = ip_graph ;
   		series x = time y = survdiff /legendlabel=' A=1 - A=0 ' lineattrs=(thickness=2) ;      
   		series x = time y = diff_2_5 /legendlabel='2.5 percentile' lineattrs=(thickness=2);  
   		series x = time  y = diff_97_5 / legendlabel='97.5 percentile' lineattrs=(thickness=2);
   		title 'Survival difference from IP weighted hazards model';
   		yaxis label = 'Survival difference'  values = (-0.1 to .1 by 0.05)  ;
   		xaxis label = 'Months '  values = (0 to 120 by 12) ;
   		keylegend/noborder ;
	run;
	title ;

	data sample0_gf ;
		set gf_graphs_boot (where = (sample = 0));
	run;

	proc sort data = gf_graphs_boot ; by time ; run;

	proc univariate data= gf_graphs_boot(where = (sample > 0)) noprint;
		var surv_diff ;
		by time ;
		output out=mysurv_pct_gf pctlpre=diff_   pctlpts=2.5,97.5; 
	run; 

	data gf_graph ;
		merge sample0_gf mysurv_pct_gf ;
		by time ;
	run;

	proc print data= gf_graph; id time; run;

	proc sgplot data = gf_graph ;
	   series x = time y= surv_diff /legendlabel=' A=1 - A=0 ' lineattrs=(thickness=2) ;      
	   series x = time y = diff_2_5 /legendlabel='2.5 percentile' lineattrs=(thickness=2);  
	   series x = time  y = diff_97_5 / legendlabel='97.5 percentile' lineattrs=(thickness=2);
	   title 'Survival difference from g-formula';
	   yaxis label = 'Survival difference'  values = (-0.1 to .1 by 0.05)  ;
	   xaxis label = 'Months '  values = (0 to 120 by 12) ;
	   keylegend/noborder ;
	run;
	title ;

	proc datasets library = work nolist ;
		delete justids_boot sample0 sample0_gf ip_graphs_boot gf_graphs_boot mysurv_gf mysurv_pct mysurv_pct_gf justids for_ipw_boot;
	quit;

%mend ;

options mprint notes nospool ;
options nomprint nonotes ;
%bootstrap(nboot = 500);





