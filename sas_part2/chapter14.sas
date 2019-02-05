
/***************************************************************
PROGRAM 14.1
Preprocessing, ranks of extreme observations, IP weights for censoring
Data from NHEFS
***************************************************************/;
libname causinf "C:\dropbox\ci\data";

/* some preprocessing of the data */
data nhefs;
    set causinf.nhefs;
    cens= (wt82 eq .);
run;

/* ranking of extreme observations */
proc univariate data= nhefs;
    id seqn;
    var wt82_71;
run;


/* estimation of denominator of ip weights for C */
proc logistic data=nhefs;
    ods exclude ClassLevelInfo Type3 Association FitStatistics GlobalTests;
    class exercise active education;
    model cens = qsmk sex race age age*age education
                smokeintensity smokeintensity*smokeintensity smokeyrs smokeyrs*smokeyrs
                exercise active wt71 wt71*wt71;
    output out=est_prob_d_c (keep= seqn pd_cens) p=pd_cens;
run;
proc sort; by seqn; run;

data nhefs_w;
    merge nhefs est_prob_d_c;
    by seqn;
    if cens=0; * observations with cens=1 only contribute to censoring models;
    w_c= 1 / pd_cens;
run;

proc univariate data=nhefs_w;
    var w_c;
    id seqn;
run;


/*****************************************************************
PROGRAM 14.2
G-estimation of a 1-parameter structural nested mean model
Brute force search
Data from NHEFS
******************************************************************/;


/*********************************************************/
/* G-estimation: Checking one possible value of psi      */
/*********************************************************/
data nhefs_gest;
    set nhefs_w;
    psi= 3.446;
    Hpsi= wt82_71-psi*qsmk;
run;

run;

proc genmod data= nhefs_gest descending;
    ods exclude ClassLevels ParmInfo ModelInfo;  
    class seqn exercise active education;
    weight w_c;
    model qsmk = sex race age age*age education
                smokeintensity smokeintensity*smokeintensity smokeyrs smokeyrs*smokeyrs
                exercise active wt71 wt71*wt71
                Hpsi/ dist= bin link= logit;
    repeated subject=seqn / type=ind;
run;
quit;



/*********************************************************/
/* G-estimation: Checking multiple possible values of psi*/
/*********************************************************/

%macro gest_snm;

    %do counter=200 %to 500 %by 10; /* integers only in do loops */
        %let psi= %eval(&counter)/100;
        data nhefs_gest;
            set nhefs_w;
            Hpsi= wt82_71 - &psi*qsmk;
        run;

        title "Testing PSI=" &psi;

        proc genmod data= nhefs_gest descending; 
            ods exclude ClassLevels ParmInfo ModelInfo;  
            class seqn exercise active education;
            weight w_c;
            model qsmk = sex race age age*age education
                        smokeintensity smokeintensity*smokeintensity smokeyrs smokeyrs*smokeyrs
                        exercise active wt71 wt71*wt71
                        Hpsi/ dist= bin link= logit;
            repeated subject=seqn / type=ind;
        run;
        quit;

    %end;

%mend;

%gest_snm;
title;


/***************************************************************
PROGRAM 14.3
G-estimation for 1-parameter structural nested mean model
Closed form estimator
Data from NHEFS
***************************************************************/;


/*************************************************************/
/* G-estimation: Closed form estimator linear mean models    */
/*************************************************************/

proc logistic data = nhefs_w descending ;
  ods exclude ClassLevelInfo Type3 Association FitStatistics GlobalTests;
  class exercise active education;
  weight w_c;
  model qsmk = sex race age age*age education
                smokeintensity smokeintensity*smokeintensity smokeyrs smokeyrs*smokeyrs
                exercise active wt71 wt71*wt71 ;
  output out = predA p = pqsmk ;

run;


proc iml ;
	use predA ;
	read all var {qsmk pqsmk wt82_71 w_c };
	close pqsmk ;
/* solve sum(w_c * H(psi) * (qsmk - E[qsmk | L]))  = 0 */
/* for a single psi and H(psi) = wt82_71 - psi * qsmk */
/* this can be solved as psi = sum( w_c * wt82_71 * (qsmk - pqsmk)) / sum(w_c * qsmk * (qsmk - pqsmk)) */
	diff = qsmk - pqsmk ;
	part1 = w_c # wt82_71 # diff ;
	part2 = w_c # qsmk # diff ;
	psi = sum(part1 ) / sum(part2) ;
	create onepsi from psi ;
	append from psi ;
	close onepsi ;
	print psi ;
quit;


 
/*************************************************************/
/* G-estimation: Closed form estimator for 2-parameter model */
/*************************************************************/

proc iml ;
    use predA ;
    read all var {qsmk pqsmk wt82_71 w_c  smokeintensity};
    close predA ;
    diff = qsmk - pqsmk ;
    diff2 = w_c # diff ;

    lhs = J(2,2,0) ;
    lhs[1,1] = sum( qsmk # diff2) ;
    lhs[1,2] = sum( qsmk # smokeintensity  # diff2 ) ;
    lhs[2,1] = sum( qsmk # smokeintensity # diff2) ;
    lhs[2,2] = sum( qsmk # smokeintensity # smokeintensity # diff2 ) ;

    rhs = J(2,1,0) ;
    rhs[1] = sum(wt82_71 # diff2 ) ;
    rhs[2] = sum(wt82_71 # smokeintensity # diff2 ) ;

    psi = t(solve(lhs,rhs));
  
  *  print lhs, rhs , psi ; 
    create twopsi from psi ;
    append from psi ;
    close twopsi ;
    print psi ;
quit;




