
/*  Roger Logan rwlogan@hsph.harvard.edu */


%macro INITIATORS(datain = ,           /* input data set containing at minimum id, period and treatment variables */
                    id = id,
                    period = period, 
                                
                    treatment = , 
                    
                    outcome = ,       

                    first_period = ,    /* first period to start expanding about */
                    last_period = ,     /* last period to expand about */

                    eligible =  ,       /* indicator of whether or not a treatment observation is eligible to be 
                                            expanded about */
                    model_switchn = ,
                    cov_switchn = ,      /* covariates to be used in logistic model for switching probabilities */
                    class_switchn = ,     /* class variables used in logistic model */

                    model_switchd = ,
                    cov_switchd = ,      /* covariates to be used in logistic model for switching probabilities */
                    class_switchd = ,     /* class variables used in logistic model */
                    
                    cense        = ,      /* censoring variable */
                    pool_cense   = 0,    /* pool the numerator and denominator models.
                                             0 = split models by previous treatment Am1 = 0 and Am1 = 1 as in treatment models
                                             1 = pool all observations together into a single numerator and denominator model
                                          */
                    model_censen = ,
                    cov_censen   = ,      /* covariates to be used in logistic model for censoring weights */
                    class_censen = ,     /* class variables used in censoring logistic model */

                    model_censed = ,
                    cov_censed   = ,      /* covariates to be used in logistic model for censoring weights */
                    class_censed = ,     /* class variables used in censorint logistic model */

                

                    lag_p_nosw = 1,      /* when 1 this will set the first weight to be 1 and use 
                                            p_nosw_d and p_nosw_n at followup-time (t-1) for calculating the 
                                            weights at followup-time t  */
                  
                                                  
                 
                    include_expansion_time = 1, /* include for_period and for_period2 in switch models */
                    include_expansion_time_case = 1, /* include for period in model for outcome */
                    include_followup_time = 1 , /* include follow up time in switch model (based on for_period ) */
                    include_followup_time_case = 1 , /* include follow up time in model for outcome */

                    include_regime_length = 0 ,                 
                    eligible_wts_0= , /* eligibility criteria used in weights for model condition Am1 = 0 */
                    eligible_wts_1= , /* eligibility criteria used in weights for model condition Am1 = 1 */

                    outcomeCov = ,     /* function of baseline covariates used in final model.   */
                    outcomeClass = ,   /* Any categorical variables used in the final model */
                    outcomeCov_var = , /* list of individual baseline variables used in final model */
                    model_var      = , /************************************************************************** 
                                           variables of interest to be used in final model. These variables need to be defined in the sub
                                           macro %create_final_variables which will use the submacro %building_blocks. These two macros are defined
                                           defined for the following two cases when model_var is left to be missing (by defualt)
                                          
                                         
                           
                                         ***************************************************************************/
                        
                    use_censor = 0,
                    first_followup = ,
                    last_followup = ,
 
                    /* analysis options */
                   
                     final_analysis = 1, /* run weighted pooled logistic regression on the person-trial data */
                     use_weights = 1  ,  /* use weights in analysis. If 0 then no weights will be calculated */
                     calculate_var = 1,  /* calculate robust variance for coefficients of final analysis models */

                     /* which types of weighted models to run when final_analysis = 1  */

                     run_unweighted_analysis = 0, /* run the final model with no weights when use_weights = 1 */
                     run_weighted_analysis = 1, /* run the final model with original weights */
                     run_p99_analysis = 0 ,     /* run the final model with truncating the weights at the 1st and 99th percentile */
                     run_user_limits_analysis = 0, /* run the final model with truncating the weights using user defined limits */
                  
                    lower_weight = ,   /* use lower weight as minimum possible weight */
                    upper_weight = ,   /* use upper weight as maximum possible weight*/
                                       
                    where_case = ,  /* where condition used in subsetting the data used in final analysis, can be multiple 
                                          conditions separated by a ":" to produce different final models */
                    where_var = , /* variables used in where_case = , need to be used when the variables are not included 
                                    in the final model */

                    run_base_model = 0 ,/* when there are where_case conditions should we run the model with no conditions */  
                    sorted = 0,      /* data sorted by id and period */
                    print_option =  , /* printing option for proc logistic either empty or noprint  */ 
                
                
                    check_missing = 0 , /* check for missing values in final model, should only be _case when use_censor = 1 */ 
                    baseline_offset = 0

                 
     );

     /* create local macro variables */
     %local by_statement missing_indicator periods_dropped rowcount ngroups block change_pooling;
     %local n_pooling_values n_pooling_times pooling_diff for_npooled  start_pooling will_change_pooling pooling_time_index ;
     %local block_name timeVarWeight controls_created cases_created;
     %local  _for_model  _n_model   for_dbeta  group for_group class_var nclass ;
     %local freq_levels num_cases ssize i0p99level i1p99level p99level i0p01level i1p01level followup_var ; 
     %local added_baseline_f added_class_f numids ;
     %let rowcount = 0 ;
     
     %let keep_list =  for_period &id  &period  switch &cov_switchd &class_switchd &cense  &cov_censed &class_censed  _assigned_treatment  &outcome _time_of_event ;
     
     %let assigned_treat = &treatment ;

   
         
      
     %if %bquote(&model_var) = %then %do;
        %if &use_censor = 0 %then %do;
            %let model_var = _dose_ _dose2_ ;            
        %end;
        %else %if &use_censor = 1 %then %do;
            %let model_var = &assigned_treat ;           
        %end;
     %end;
      

     %if %eval(&include_expansion_time) = 1 %then %let keep_list = &keep_list for_period2 ;
     %if %eval(&include_followup_time)  = 1 %then %let keep_list = &keep_list _followup_time _followup_time2;
     %let for_group = _group ;
   
     %let keep_time_var = &id for_period for_period2 &period &assigned_treat &outcomeCov_var &outcomeClass  &for_group  ;
     %if %bquote(&eligible)= %then %let eligible = _eligible ;

     %if %eval(&use_weights) = 0 %then %do;
            %let run_weighted_analysis = 0;
            %let run_p99_analysis = 0 ;
            %let run_user_limits_analysis = 0 ;
           
     %end;
     
   %if %bquote(&where_case)^= %then %do;
      %let case_cond = 0 ;
      %let _holder_ = &where_case ;
      %let n = 1;
      %do %until(%qscan(&&_holder_,%eval(&n),%str(:)) = %str());
      %let where_cond&n = %qscan(&_holder_,%eval(&n),%str(:));
      %let n = %eval(&n + 1) ;
      %end;
      %let case_cond = %eval(&n - 1);

      %put case_cond = &case_cond ;
      %do ii = 1 %to %eval(&case_cond);
      %put condition &ii is &&where_cond&ii ;
      %end;
      %let keep_list      = %remove_duplicate(&keep_list &where_var) ;
      %let keep_time_var  = %remove_duplicate(&keep_time_var &where_var) ;  
      %let outcomeCov_var = %remove_duplicate(&outcomeCov_var &where_var) ;
   %end;
   %else %do;
      %let case_cond = 0;
      %let where_var = ;      
      %if %eval(&final_analysis) = 1 %then %let run_base_model = 1;
      
   %end;
   
   
    %if &run_user_limits_analysis = 0 %then %do;
         %let lower_weight = ;
         %let upper_weight = ;
    %end;   
    
    %let followup_var = _followup_time _followup_time2 ;
   
    %if %eval(&use_weights) = 0 & %eval(&final_analysis) = 1 %then %let run_unweighted_analysis = 1 ;
 
    %initialize_datasets ;


 

    proc means data = _datain_ (keep =  &period  )  min max  noprint;
    var  &period ;         
    output out = periodrange max = maxperiod min = minperiod;
    run;

    data _null_;
    set periodrange;
    call symput('period_max',maxperiod);
    call symput('period_min',minperiod);
    run;


    %if %bquote(&first_period) = %then %let first_period = %eval(&period_min);                 
    %if %bquote(&last_period)  = %then %let last_period  = %eval(&period_max);
        
    %if %eval(&first_period) > %eval(&period_min) %then %do;
        data work._datain_;
        set work._datain_;
        if &period >= %eval(&first_period) ;
        run; 
    %end;

    %let  period_num = %eval(&period_max - &period_min + 1);


    /* run main analysis */

    %method5 ;

      
    proc datasets library=work   nolist;
    delete _datain_  Periodrange                      ;
    quit;
      
%mend ;

%macro initialize_datasets ;
        %let ncov = %numargs(&outcomeCov_var);
        data _datain_  %if &ncov > 0 %then (drop = _length_) ;;
        set &datain (keep = &id &period  &treatment &outcome &eligible_wts_0 &eligible_wts_1 
                            &cov_switchd &class_switchd
                           %if &eligible ^= _eligible_  %then &eligible ;
                           %if %bquote(&cense)^= %then &cense  &cov_censed &class_censed ;                         
                           %if %eval(&final_analysis) = 1   %then &outcomeCov_var; /* &where_var ; */                       
                    )  ;

        %if &eligible = _eligible_ %then _eligible = 1 ;;
        if _n_ = 1 then do ;
          %do i = 1 %to &ncov ;
               %let word0 = %scan(&outcomeCov_var,%eval(&i),%str( ));
               _length_ = vlength(&word0);
               %global length&i ;
               call symput("length&i",trim(left(_length_)));
          %end;
       end; 
       run;

       %if %eval(&sorted ) = 0 %then %do;

          proc sort data = work._datain_ ;
          by &id &period;
          run;
       %end;

  

      data _for_outcome_;
      set _datain_;
      length _time_of_event 4 ;
      by &id &period;
      if last.&id;
      _time_of_event = 9999;
      if &outcome = 1 then _time_of_event = &period ;
      keep &id _time_of_event ;
      run;

      data _datain_;
      merge _datain_    _for_outcome_;
      by &id;
      run;
         

      %let class_var =  %nrstr(&outcomeClass);
      %let nclass = 1 ;
      %let cov_var = %nrstr( &outcomeCov) ;
      %let aaaa = %unquote(&class_var) ;
      %if %bquote(&aaaa) ^=  %then  %remove_extraneous;
       
      data _datain_ (drop = &id );
      set _datain_ end = _end_;
      retain  _regime_start  _newid_ _cumA_ ; 
      by &id ;
      if _n_ = 1 then _newid_ = 0 ;
      length switch   3 ;
                           
      _Am1_ = lag(&treatment);
              
       if first.&id then do;
             _newid_ = _newid_ + 1;
             _cumA_ = 0;
             _Am1_ = 0 ;
             switch = 0 ; /* set there to be no switch at "baseline" observation */
             _regime_start = &period ;
             _time_on_regime_ = 0;
             _time_on_regime2_ = 0 ;                     
       end;
       else do ;
            if _Am1_ ^= &treatment then switch = 1;
            else switch = 0 ;                 
            _time_on_regime_ = &period - _regime_start  ;
            _time_on_regime2_ = _time_on_regime_ ** 2 ;
            if switch = 1 then _regime_start = &period ;
       end;

       _cumA_ = _cumA_ + &treatment ;  /* for dose variable in method 5 */                        
       if _end_ then call symput("numids",trim(left(_newid_)));
       run;   

       %let id = _newid_ ;         
  
   

%mend;


%macro numargs(arg);

     %let n = 1;
     %if &arg^= %then %do;
          %do %until (%scan(&arg,%eval(&n),%str( ))=%str());
               %let word = %scan(&arg,&n);
               %let n = %eval(&n+1);
          %end;
     %end;
     %eval(&n-1) /* there is no ; here since it will be used as %let a = %numargs(&b) ;
     and the ; is included at the end of this line  */                            
                            
%mend numargs;


   
%macro robust(outest=a, out=_last_, id=id, df=);

     /*********************************************************************
     **********************************************************************
     MACRO ROBUST produces robust estimates of the covariance matrix for
     the coefficients from PROC LOGISTIC or PROC PHREG when the data are
     clustered. One common application is to longitidunal data with each
     individual treated as a cluster.  IML is required.

     The macro uses the method of Halbert White (1982) "Maximum likelihood
     estimation of misspecified models," Econometrica 50: 1-25.

     Author:  Paul D. Allison, University of Pennsylvania
     allison@ssc.upenn.edu

     Adapted from an IML program written by Terry Therneau, Mayo Clinic.

     For PROC LOGISTIC, you must specify the OUTEST=name1 and
     COVOUT options on the PROC statement. You must also use the OUTPUT
     statement with OUT=name2 and DFBETAS=namelist.  The namelist should
     have one name for each term in the model, including the intercept.
     There must also be a variable in the data set containing a unique
     value (either character or numeric) for each cluster. The macro is
     invoked after running the LOGISTIC procedure.

     For PROC PHREG, you must specify OUTEST=name1 on the PROC
     statement. (COVOUT is unncessary).  You must also use the OUTPUT
     statement with OUT=name2 and DFBETA=namelist.  The namelist should
     have one name for each variable in the model. There must also be a
     variable in the data set containing a unique value (either character
     or numeric) for each cluster. This variable must be added to the
     OUTPUT data set by using the ID statement.  The macro is
     invoked after running the PHREG procedure.

     The macro has the following parameters:

     OUTEST   Name of data set used in the OUTEST= option.
     OUT      Name of data set used in the OUT= option.
     ID       Name of variable containing a unique value for each cluster
     DF       List of names used in the DFBETAS or DFBETA option.

     Examples of usage:

     proc logistic outest=a covout;
     model y = x z w;
     output out=b dfbetas=dint dx dz dw;
     run;
     %robust(outest=a, out=b, id=subjid, df=dint dx dz dw)

     proc phreg outest=a;
     model y*d(0) = x z w;
     id subjid;
     output out=b dfbeta= dx dz dw;
     run;
     %robust(outest=a, out=b, id=subjid, df=dx dz dw)

     BE CAREFUL: it's DFBETAS in LOGISTIC but DFBETA in PHREG.
     (The former is standardized, the latter is not).
     Also PHREG does NOT have an intercept.
     **************************************************************
     **************************************************************/


     proc means data=&out noprint;
     class &id;
     var &df;
     output out=_out1_(keep=&df) sum=&df;
     types &id ;
     run;


     data _out1_;
     set _out1_ ;
     array d (*) &df;
     if d(1)=. then delete;
     do _i_ = 1 to dim(d) ;
         if d(_i_) = . then d(_i_) = 0 ;
     end;
     drop _i_ ;
     run;

     data _reduce_;
     set &outest;
     array abc(*) _character_;
     length name $8;
     call vname(abc(1),name);
     if name ne '_LINK_' and _type_ eq 'COV' then delete;
     drop  _lnlike_;
     run;


     proc iml;
     use _reduce_ where (_type_='COV');
     read all into cov;
     use _reduce_ where (_type_='PARMS');
     read all into b[colname=vname];
     if ncol(cov)=0 then se=1;
     else se=sqrt(diag(cov));

     nn = ncol(se);
     do i = 1 to nn ;
        if se[i,i] = . then se[i,i] = 0 ;
     end;

     use _out1_;
     read all into x;
     x=x*se;
     v=x`*x;
     se=sqrt(vecdiag(v));
     reset noname fuzz=.000001;

     se = t(se);
     create se_block from se ;
     append from se ;

     create covar_block from v;
     append from v ;

     quit;
     run;
   

%mend robust;
 
%macro remove_extraneous;
  
     %do i = 1 %to &nclass;

          %let _var = %qscan(&class_var, &i,%str( ));
          %let _cvar = %qscan(&cov_var,&i,%str( ));
          %let _class = %unquote(&_var);
          %let _allcov = %unquote(&_cvar);
          %if %bquote(&_class)^= %then %do; 

               /* remove entries in list1 that are not in list2 */
            %local   test word ;
            %let ii = 1;
            %let _class_used= ;
            %do %until(%scan(&_class,&ii,%str( )) = %str() );

                %let test = %scan(&_class,&ii);
                %let j = 1;
                %let remove = 1;
                %do %until(%scan(&_allcov,&j,%str( ))=%str());
                     %let word = %scan(&_allcov,&j);
                      %if &test = &word %then %let remove = 0 ;
                      %let j = %eval(&j + 1);
                %end;
                %if %eval(&remove) = 0 %then %let _class_used = &_class_used &test ;;
                %let ii = %eval(&ii + 1);

            %end;

           %let %substr(&_var,2) = &_class_used;

        %end;  
             
     %end;

%mend remove_extraneous;

%macro remove_duplicate(list1) ;
     %local leng1 i word word2 j list_tmp good ;
       %let leng1=%numargs(&list1);

       %let list_tmp = %scan(&list1,1,%str( ));
       %do i = 2 %to &leng1 ;
           %let word = %scan(&list1,%eval(&i),%str( ));
           %let good = 1 ;
           %do j = 1 %to %eval(&i-1);
                %let word2 = %scan(&list_tmp,%eval(&j),%str( ));
                %if %cmpres(&word2)=%cmpres(&word) %then %let good = 0 ;
           %end;
           %if &good = 1 %then %let list_tmp = &list_tmp &word ;
       %end;
       &list_tmp 
%mend ;


 
%MACRO RCSPLINE(x,knot1,knot2,knot3,knot4,knot5,knot6,knot7,
                  knot8,knot9,knot10, norm=2);

/* MACRO RCSPLINE

   For a given variable named X and from 3-10 knot locations,
   generates SAS assignment statements to compute k-2 components
   of cubic spline function restricted to be linear before the
   first knot and after the last knot, where k is the number of
   knots given.  These component variables are named c1, c2, ...
   ck-2, where c is the first 7 letters of X.

   Usage:

   DATA; ....
   %RCSPLINE(x,knot1,knot2,...,norm=)   e.g. %RCSPLINE(x,-1.4,0,2,8)

        norm=0 : no normalization of constructed variables
        norm=1 : divide by cube of difference in last 2 knots
                 makes all variables unitless
        norm=2 : (default) divide by square of difference in outer knots
                 makes all variables in original units of x

   Reference:

   Devlin TF, Weeks BJ (1986): Spline functions for logistic regression
   modeling. Proc Eleventh Annual SAS Users Group International.
   Cary NC: SAS Institute, Inc., pp. 646-51.


   Author  : Frank E. Harrell Jr.
             Clinical Biostatistics, Duke University Medical Center
   Date    : 10 Apr 88
   Mod     : 22 Feb 91 - normalized as in S function rcspline.eval
             06 May 91 - added norm, with default= 22 Feb 91
             10 May 91 - fixed bug re precedence of <>

                                                                      */


   %LOCAL j v7 k tk tk1 t k1 k2;
   %LET v7=&x; %IF %LENGTH(&v7)=8 %THEN %LET v7=%SUBSTR(&v7,1,7);
     %*Get no. knots, last knot, next to last knot;
   %DO k=1 %TO 10;
       %IF %QUOTE(&&knot&k)=  %THEN %GOTO nomorek;
   %END;
   %LET k=11;
   %nomorek: %LET k=%EVAL(&k-1); %LET k1=%EVAL(&k-1); %LET k2=%EVAL(&k-2);
   %IF &k<3 %THEN %PUT ERROR: <3 KNOTS GIVEN.  NO SPLINE VARIABLES CREATED.;
   %ELSE %DO;
      %LET tk=&&knot&k;
      %LET tk1=&&knot&k1;
      DROP _kd_; _kd_=
      %IF &norm=0 %THEN 1;
      %ELSE %IF &norm=1 %THEN &tk - &tk1;
      %ELSE (&tk - &knot1)**.666666666666; ;
      %DO j=1 %TO &k2;
         %LET t=&&knot&j;
         &v7&j=max((&x-&t)/_kd_,0)**3+((&tk1-&t)*max((&x-&tk)/_kd_,0)**3
             -(&tk-&t)*max((&x-&tk1)/_kd_,0)**3)/(&tk-&tk1)%STR(;);
       %END;
   %END;
%MEND;



%macro method5;

     %local model_suffix first_model number_of_models number_of_passes first_model ;

     %let class_holder = &outcomeClass;
     %let cov_holder = &outcomeCov_var;
    
     %if %eval(&case_cond) > 0 %then %do;
         %let number_of_passes = %eval(&case_cond) ;
      %if %eval(&run_base_model) = 1 %then %let number_of_passes = %eval(&number_of_passes + 1);
     %end;
     %else %let number_of_passes = 1 ;

     %do pass = 1 %to &number_of_passes ;
     %local i0p01level&pass i1p01level&pass p01level&pass i0p99level&pass i1p99level&pass p99level&pass;
     %end;
      

     %let number_of_models = %eval(&run_unweighted_analysis + &run_weighted_analysis + &run_p99_analysis + &run_user_limits_analysis);
    
     %method5_analysis ;
        
     %if  &calculate_var = 1  %then %do;
          
        /***
         unweighted analysis   ==> model = 0
         weighted analysis     ==> model = 1
         p99 weights           ==> model = 2
         user limited weights  ==> model = 3

         ****/
     
         %let all_models = &run_unweighted_analysis &run_weighted_analysis &run_p99_analysis &run_user_limits_analysis ;


         %if &use_weights = 0 %then %let number_of_models = 0 ; /* the loop below will start at 0 and not 1 */

         %if %eval(&final_analysis) = 1   %then %do;
        
              %do pass = 1 %to %eval(&number_of_passes) ;


                  %do model = 0 %to 3;
                      /* if only running unweighted analysis and a truncated weight model, the sequence of
                         model values is 0, 2, and 3 and number of models = 3 . need to skip model = 1 */

                      %let run_model = %scan(&all_models,%eval(&model + 1),%str( ));

                      %if &run_model = 1 %then %do;
                   
               
                     %let suffix = &pass._&model ;
                               
                           proc iml;
                           use _se&suffix ;
                           read all into se;
                           use _est_betas&suffix;
                           read all into betas;
                           avg_var = betas // se ;
                           create _avg from avg_var;
                           append from avg_var ;
                           quit;
 
                           proc transpose data = _avg out = output&suffix (drop = _NAME_);
                           run;
   
                           data _tmp;
                           set _est_betas&suffix (obs = 1 );
                           run;

                           proc transpose data = _tmp out= _tmpt (drop = COL1);
                           run;

                           data _tmpt ;
                           set _tmpt;
                           if _LABEL_ = '' then _LABEL_ = _NAME_;
                           if _n_ = 1 then _LABEL_ = 'intercept';
                           run;
        
                          data output&suffix ;
                          merge output&suffix (rename = (   COL1 = estimate COL2 = std)) _tmpt (drop = _NAME_ );
                          rename _LABEL_ = name ;
  
                          label p_value = 'Pr > |Z|';
                          format p_value PVALUE8.4 ;  
                          if std > 0 then do ;   
                              lb = estimate - 1.96 * std;
                              ub = estimate + 1.96 *std ;
                              z = estimate / std;
                              p_value = 2*(1-probnorm(abs(z))); 
                          end;
                          else std = . ;
                          run;
            


                      %if %eval(&pass) = 1 & %eval(&run_base_model) = 1 %then %let where_case_pass = ;
                         %else %do;
                             %let _tmp = %eval(&pass - 1);
                             %if &run_base_model = 0 %then %let where_case_pass = &&where_cond&pass;
                             %else %if &run_base_model=1 %then  %let where_case_pass = &&where_cond&_tmp; 
                         %end;
                         %if %eval(&model) = 0 %then %do;
                       %if %bquote(&where_case_pass) = %then %do;             
                                  title "Analysis with no weights and using robust variance";
                       %end;
                       %else %do;
                           title "Analysis with no weights where &where_case_pass using robust variance";
                             %end;      
                         %end;
                         %else %if %eval(&model) = 1 %then %do;
                        %if %bquote(&where_case_pass) = %then %do;            
                                   title "Analysis  with original weights and using robust variance";
                        %end;
                        %else %do;
                             title "Analysis with original weights where &where_case_pass using robust variance";
                              %end;          
                        %end;                 
                        %else  %do ;
                             %if %bquote(&where_case_pass) = %then %do;                       
                                  title1 "Analysis using robust variance and truncating weights to be between  ";
                       %end;
                       %else %do;
                            title1 "Analysis using robust variance when &where_case_pass and  truncating weights to be between ";
                       %end;                              
                             %if &model = 2 %then %do;
                                  title2 "&&p01level&pass (first percentile) and &&p99level&pass (99th percentile). ";               
                             %end;
                             %else %if &model = 3 %then %do;
                                  title2 "user defined limits of &lower_weight and &upper_weight ";
                             %end;  

                       %end;


                       proc print data = output&suffix noobs;
                       var name estimate std lb ub z p_value ;
                       run;
                       title ;

         
                       proc sql noprint;
                       select _NAME_
                       into :_for_model separated by ' '
                       from  _tmpt;
                       quit;

                       %let _n_model = %numargs(&_for_model) ;

         
                      data _se&suffix ;
                      set _se&suffix;
                      rename  
                         %do i = 1 %to &_n_model;
                            %let vartmp = %scan(&_for_model,%eval(&i),%str( ));                
                            COL&i = &vartmp 
                         %end;
                      ;
                      run;

                     data _covar&suffix ;
                     set _covar&suffix ;
                     rename  
                         %do i = 1 %to &_n_model;
                            %let vartmp = %scan(&_for_model,%eval(&i),%str( ));                
                            COL&i = &vartmp 
                         %end;
                      ;
                      length name $20 ;
                   
                      %do i = 1 %to &_n_model;
                            %let vartmp = %scan(&_for_model,%eval(&i),%str( ));                
                            if _n_ = &i then name  = "&vartmp" ;
                      %end;
                
                     run;
 
                     %if &use_censor = 0 %then %do ;
                         proc print data = _covar&suffix (where= (index("&model_var",compress(name)) > 0)    );
                         var name &model_var;
                         run;
                    %end; 
                  %end ; /* run_model = 1 */
              %end; /* model */
            
            %end;
         %end; /* pass */
   %end;

%mend;    


%macro method5_analysis;

    %local models _num_models  first_model;

    /* calculate weights and create person-trial dataset for final analysis */

    %cw3   ;

/* need to fix these two variables */
   %let group = 1;
   %let first_group = 1 ;
 
  

     /** order of possible models 
          1) pooled logistic model with no weights
                   further break down depending on number of different models in where_case variable 
          2) pooled logistic model using weights 
               further break down depending on the number of conditions given in the where_case macro variable
          3) pooled logistic model using weights with maximum being set by the 99th percentile
               
          4) additional model based on values of lower_weight and lower_weight, when p99 = 1  

     **/ 

    

    %if %eval(&final_analysis) = 1   %then %do;

      
       %if %bquote(&first_followup) = %then %let first_followup = 0;
       %if %bquote(&last_followup) =  %then %let last_followup = %eval(&period_max + 1);       
   

       %let weight_list = _weight_;
       %if &run_p99_analysis = 1 %then %let weight_list = &weight_list _weight_p99;
       %if %bquote(&lower_weight)^= %then %let weight_list = &weight_list _weight_pX ;
       
   
      /* pass variable indicates number of subset analyses to run for each model */

      %do pass = 1 %to %eval(&number_of_passes) ;

      
      %if %eval(&pass) = 1 & %eval(&run_base_model) = 1 %then %let where_case_pass = ;
      %else %do ;     
              %let _tmp = %eval(&pass - 1);  
              %if %eval(&case_cond) = 0 %then %let where_case_pass =  ;             
              %else %if &run_base_model = 1 %then %let where_case_pass = &&where_cond&_tmp;
              %else %if &run_base_model = 0 %then %let where_case_pass = &&where_cond&pass ; 
       %end;
     
                   

      %if %bquote(&where_case_pass)^= %then title1 "Weights for those subjects with &where_case_pass";;
      title2 "Analysis of weights for switching treatment by pooling &assigned_treat values";
          



         %if &run_p99_analysis = 1  or &run_user_limits_analysis = 1 %then %do;
         
             %if &run_p99_analysis = 1  %then %do;
                 proc univariate data = _switch (keep =&id  _weight_ _assigned_treatment &where_var  rename = (_assigned_treatment = &assigned_treat));
                 id &id ;
                 %if %bquote(&where_case_pass)^= %then where &where_case_pass ;;                
                 var _weight_ ;              
                 output out = _p99_ p99 = p99  p1 = p1 ;
                 run;
        
                  
                 data _p99_;
                 set _p99_;                  
                 call symput("p99level&pass" , p99);
                 call symput("p01level&pass" ,  p1);
                 run;
             %end;

          data _switch_pass ;
          set _switch %if %bquote(&where_case_pass) ^= %then (where = ( &where_case_pass )) ;;
             %if &run_p99_analysis = 1 %then %do;
               
                  _weight_p99 = _weight_ ;   
                
                     if _weight_ > &&p99level&pass then _weight_p99 = &&p99level&pass ;
                     if _weight_ < &&p01level&pass then _weight_p99 = &&p01level&pass ; 
             %end;  
     
          %if &run_user_limits_analysis = 1 %then %do;
               _weight_pX = _weight_ ;                   
                     if _weight_ > &upper_weight then _weight_pX = &upper_weight ;
                     if _weight_ < &lower_weight then _weight_pX = &lower_weight ;
          %end;
             run ;

         %end ;
         %else %do;

             proc means data = _switch(keep =  _weight_ _assigned_treatment rename = (_assigned_treatment = &assigned_treat))
                     n min p1 p10 p25 p50 p75 p90 p99 max mean var std  ;
             var _weight_ ;
             run;

            data _switch_pass ;
            set _switch %if %bquote(&where_case_pass) ^= %then (where = ( &where_case_pass )) ;;
            run;



            
         %end;


        %if &check_missing=1 %then %do;
            /*********** 
            check for missing values in variables in final models, when use_censor = 1 only missing
            values should be in _case variable due to censoring at a change in treatment. When 
            use_censor = 0 there should be no missing values ( what about original missing outcome
            due to censoring in original data
             ***********/
            proc means data = _switch_pass n nmiss ;
            var  &outcome  &model_var &outcomeCov  
                      %if &include_expansion_time_case = 1 %then for_period for_period2 ;
                      %if %eval(&include_followup_time_case) = 1 %then &followup_var ; ; 
           run;
        %end;

        %if &run_unweighted_analysis = 1 %then %do;   
             
              %let model = 0 ;
              %let drop_list = ;     
              %let drop_list =  &outcome   &model_var  for_period for_period2  &followup_var &outcomeCov_var  ;
                      
               
              proc logistic 
                  data = _switch_pass ( keep = &id  &outcome &model_var &outcomeCov_var _followup_time  /* &where_var*/
                                                %if &include_expansion_time_case = 1 %then for_period for_period2 ;
                                                %if %eval(&include_followup_time_case) = 1 %then &followup_var ;                                           
                                         where = ( %eval(&first_followup) <= _followup_time <= %eval(&last_followup)))
                          descending
                          outest = _est_block covout ;
              ods exclude ClassLevelInfo Type3 Association  /* GlobalTests */ Oddsratios;
              
              %if %bquote(&outcomeClass)^= %then class &outcomeClass   ;;
                   
              %if %bquote(&where_case_pass)^= %then %do;                       
                        title  "Analysis using  no weights where &where_case_pass";
              %end;
              %else %do;
                title "Analysis using no weights";
              %end;
          
              model  &outcome =  &model_var                      
                             %if %eval(&include_followup_time_case) = 1 %then &followup_var ; 
                             %if &include_expansion_time_case = 1 %then for_period for_period2 ;
                             &outcomeCov 
                            / nologscale  ridging=none
                             ; 

              %if &calculate_var = 1 %then %do;
                    output out =_for_robust(   drop = &drop_list )   dfbeta = _ALL_  ;
              %end;    
              run;
   
              data _est_tmp;
              set _est_block ;
              if _n_ = 1 ;                            
              drop _LINK_ _TYPE_ _STATUS_ _NAME_ _LNLIKE_ ;
              run;
          
              %if &calculate_var = 1 %then %do;                        
 
                   data _for_dfbetas;
                   set _for_robust(obs = 1);
                   run;

                   proc transpose data = _for_dfbetas out = _for_dfbetast;
                   run;

                   data _for_dfbetast;
                   set _for_dfbetast (keep = _NAME_);
                   if lowcase(_NAME_) = lowcase("&id") then delete;
                   run;

                   proc sql noprint;
                   select _NAME_
                   into :for_dbeta separated by ' '
                   from _for_dfbetast;
                   quit;

                    %robust(outest=_est_block, out=_for_robust, id=&id, df=&for_dbeta);
               
                 /* output of robust macro is a data set named se_block */


                   %if %eval(&group)= %eval(&first_group) %then %do;           
                        data _se&pass._&model ;
                        set se_block ;
                        run;

                        data _est_betas&pass._&model;
                        set _est_tmp;
                        run;


                        data _covar&pass._&model ;
                        set covar_block;
                        run;

                        proc datasets library = work nolist;
                        delete se_block   _for_robust _est_block _reduce_ _out1_  ;
                        quit;
                   %end;
                   %else %do;
                         proc datasets library = work nolist;
                         append base = _se&pass._&model        data = se_block;
                         append base = _est_betas&pass._&model data = _est_tmp ; 
                         append base = _covar&pass._&model     data = covar_block ;
                         delete se_block _est_tmp _for_robust _est_block _reduce_ _out1_ &group_name.&group covar_block ;
                         quit;
                    run;
                   %end;
        
               %end;  
       %end;
                           
       %let number_of_weighted_models = %eval(&run_weighted_analysis + &run_p99_analysis + &run_user_limits_analysis);
       %let weighted_models = &run_weighted_analysis &run_p99_analysis &run_user_limits_analysis ; /* ( X Y Z ) list */

       %do model = 1 %to 3  ;
          
                %let run_model = %scan(&weighted_models,%eval(&model),%str( ));

                %if &run_model = 1 %then %do;
                
                    %if %eval(&model) = 1 %then %do;
                  %if %bquote(&where_case_pass) = %then %do;             
                            title "Analysis  using original weights";
                  %end;
                  %else %do;
                     title "Analysis using original weights where &where_case_pass";
                        %end;      
                   %end;                 
                   %else  %do ;
                        %if %bquote(&where_case_pass) = %then %do;                       
                            title1 "Analysis when truncating weights using";
                  %end;
                  %else %do;
                      title1 "Analysis when truncating weights where &where_case_pass and using";
                  %end;
                  %if &model = 2 %then %do;
                      title2 "lower limit of &&p01level&pass and upper limit of &&p99level&pass ";
                  %end;
                  %else %if &model = 3 %then %do;
                      title2 "a user defined lower limit of &lower_weight and upper limit of &upper_weight";
                  %end;
                    
               %end;
           
               %let drop_list = ;     
               %let drop_list = &outcome    &model_var   for_period for_period2  &followup_var &outcomeCov_var ;
            
               proc logistic 
                      data = _switch_pass( keep = &id   &outcome  &model_var /* &where_var*/
                                             %if       &model = 1   %then  _weight_ ;
                                             %else %if &model = 2 %then _weight_p99 ;
                                             %else %if &model = 3 %then _weight_pX ;
                                             &outcomeCov_var  
                                             %if &include_expansion_time_case = 1 %then for_period for_period2 ;
                                             %if &include_followup_time_case = 1 %then &followup_var  ;
                                       where = ( &first_followup <= _followup_time <= &last_followup )
                                              )
                      descending  
                      outest = _est_block covout  ;
                  
                   ods exclude ClassLevelInfo Type3 Association  /* GlobalTests */ Oddsratios;
                   %if %bquote(&outcomeClass)^= %then class &outcomeClass   ;;

                    /* where statement is included in the data set construction above */
                     
                   model  &outcome  = &model_var
                                %if %eval(&include_followup_time_case) = 1 %then &followup_var  ;
                                %if &include_expansion_time_case = 1 %then for_period for_period2 ;                               
                                 &outcomeCov 
                               / nologscale ridging=none
                     
                     ;    
                   weight  %if  %eval(&model) = 1   %then  _weight_ ;
                          %else %if %eval(&model) = 2 %then _weight_p99 ;
                          %else %if %eval(&model) = 3 %then _weight_pX ;
                   
                          ; 
                   %if &calculate_var = 1 %then %do;
                    output out =_for_robust(  drop =    %if       &model = 1   %then  _weight_ ;
                                                           %else %if &model = 2 %then _weight_p99 ;
                                                           %else %if &model = 3 %then _weight_pX ; 
                                                           &drop_list ) dfbeta = _ALL_   ;
                   %end;
                   run;

                  data _est_tmp;
                  set _est_block ;
                  if _n_ = 1 ;                                
                  drop _LINK_ _TYPE_ _STATUS_ _NAME_ _LNLIKE_ ;
                  run;
            
                  %if &calculate_var = 1 %then %do;                        
 
                    data _for_dfbetas;
                    set _for_robust(obs = 1);
                    run;

                    proc transpose data = _for_dfbetas out = _for_dfbetast;
                    run;

                    data _for_dfbetast;
                    set _for_dfbetast (keep = _NAME_);
                    if lowcase(_NAME_) = lowcase("&id") then delete;
                    run;

 
                    proc sql noprint;
                    select _NAME_
                    into :for_dbeta separated by ' '
                    from _for_dfbetast;
                    quit;

                    %robust(outest=_est_block, out=_for_robust, id=&id, df=&for_dbeta);
               
                    /* output of robust macro is a data set named se_block */

                    %if %eval(&group)= %eval(&first_group) %then %do;           
                       data _se&pass._&model ;
                       set se_block ;
                       run;

                       data _est_betas&pass._&model;
                       set _est_tmp;
                       run;

                       data _covar&pass._&model ;
                       set covar_block ;
                       run;

                       proc datasets library = work nolist;
                       delete se_block   _for_robust _est_block _reduce_ _out1_ _group1 ;
                       quit;
                   %end;
                   %else %do;
                       proc datasets library = work nolist;
                       append base = _se&pass._&model        data = se_block;
                       append base = _est_betas&pass._&model data = _est_tmp ; 
                       append base = _covar&pass._&model     data = covar_block ;
                       delete se_block _est_tmp _for_robust _est_block _reduce_ _out1_ &group_name.&group covar_block;
                       quit;
                       run;
                   %end;


                  %end; /* calculate var = 1 */
                  title ;
                %end ; /* run_model = 1 */  
              %end; /* end of model loop */
     
        
          
          proc datasets library = work nolist ;
          delete _switch_pass ;
          quit; 
     %end; /* end pass loop */
  %end ;  /* final analysis = 1 */
%mend;



%macro cw3 ;
     
     data work.sw (index = (id_period = (&id &period))) ;
     set _datain_ ;
     by &id &period ;
     %if &use_censor=1 %then retain  _started0 _started1 _stop0 _stop1 ;;
     length _eligible0 _eligible1 3  %if &use_censor = 1 %then _started0 _started1 _stop0 _stop1 /* _starttime */ _eligible0_sw _eligible1_sw 3 /* _starttime2 5 */;  ;

     if first.&id = 1 then do;
         
          %if &use_censor = 1 %then %do;
             _started0 = 0 ;
             _started1 = 0 ;
             _stop0 = 0;
             _stop1 = 0 ;
             _eligible0_sw = 0 ;
             _eligible1_sw = 0 ;
           
          %end;
        
     end;
   
    _eligible0 = 0 ;
    _eligible1 = 0;

    %if &use_censor = 1 %then %do;    
         
 /* we stopped the follow-up on the previous observation. reset the variables to  check for a new follow-up */
 
           if _stop0 = 1 or _stop1 = 1 then do ;
               _stop0 = 0;
               _stop1 = 0 ;
               _started0 = 0 ;
               _started1 = 0 ;
               _eligible0_sw = 0 ;
               _eligible1_sw = 0 ;

                      
            end;

          /* not following any kind of follow-up yet */

           if _started0 = 0 & _started1 = 0  & &eligible = 1 then do ;
                if      &treatment = 0 then       _started0 = 1 ;                     
                else if &treatment = 1 then       _started1 = 1;
                           
           end;
    
           if _started0 = 1 and _stop0 = 0 then do ;
               _eligible0_sw = 1 ;
               _eligible1_sw = 0 ;
           end;
           else if _started1 = 1 and _stop1 = 0 then do ;
               _eligible0_sw = 0 ;
               _eligible1_sw = 1 ;
           end;
           else do ;
             _eligible0_sw = 0;
             _eligible1_sw = 0;
           end;
    

          /* can not switch treatment on first observation */

           if    switch = 1 then do ;
 
              /* on the current observation _am1_  and current treat are different. It is
                 possible that we can start a new trial at the current time if eligible. */
                 if &eligible = 1 then do ; 
                    /* start a new trial here */
                    if &treatment = 1 then do ;
                         _started1 = 1 ;
                         _stop1    = 0 ;
                         _started0 = 0 ;
                         _stop0    = 0 ;                       
                         _eligible1_sw = 1;
                    end;
                    else if &treatment = 0 then do ;
                         _started0 = 1;
                         _stop0    = 0;
                         _started1 = 0;
                         _stop1    = 0;
                         _eligible0_sw = 1 ;
                    end;
              end;
              else do ;
                 /* current observation is not eligible to start a new trial so the 
                    previous trials must stop here */                              
               _stop0 = _started0 ;
               _stop1 = _started1 ; 
              end;
         end; /* switch = 1 */
         
          if _eligible0_sw = 0 and _eligible1_sw = 0 then delete ;
 
    %end;
   
    if      _Am1_ = 0 then _eligible0 = 1 ;
    else if _Am1_ = 1 then _eligible1 = 1 ;
   
    run ;

    %if &include_regime_length = 1   %then %do;
          %let model_switchd = &model_switchd _time_on_regime_  _time_on_regime2_ ;
          %let model_switchn = &model_switchn _time_on_regime_  _time_on_regime2_ ;
    %end;

    %if &use_weights = 1 %then %do;

    proc logistic data =  work.sw(where = ( _eligible0 = 1 %if %bquote(&eligible_wts_0)^= %then and &eligible_wts_0 = 1 ;)) descending    &print_option ;
    ods select ModelInfo ResponseProfile ParameterEstimates ;                          
    title3 "Model for P(&treatment = 1 |  &treatment = 0 ,use_censor = &use_censor ) for denominator ";
   
    %if %bquote(&class_switchd)^= %then class &class_switchd  &added_class ;;
    model &treatment  = &model_switchd      / nologscale   ;
    %if %bquote(&by_statement)^= %then by &by_statement ;; 
    output out =  work._switch_d0    (keep =  &id &period  p0_d _eligible0  ) p = p0_d ;
    run; 

     proc logistic data =  work.sw(where = ( _eligible0 = 1  %if %bquote(&eligible_wts_0)^= %then and &eligible_wts_0 = 1 ; )) 
           descending    &print_option ;
     ods select ModelInfo ResponseProfile ParameterEstimates ; 
     title3 "Model for P(&treatment = 1 |  &treatment = 0 ,use_censor = &use_censor ) for numerator ";
     
     %if %bquote(&class_switchn)^= %then class &class_switchn &added_class ;;
     model &treatment = &model_switchn   / nologscale     ;
     %if %bquote(&by_statement)^= %then by &by_statement ;; 
     output out =  work._switch_n0   (keep =  &id &period     p0_n _eligible0  ) p = p0_n ;
     run; 

     proc logistic data =  work.sw(where = ( _eligible1 = 1 %if %bquote(&eligible_wts_1)^= %then and &eligible_wts_1 = 1 ; )) descending    &print_option ;
     ods select ModelInfo ResponseProfile ParameterEstimates ; 
     title3 "Model for P(&treatment = 1 |  &treatment = 1 ,use_censor = &use_censor ) for denominator ";
   
     %if %bquote(&class_switchd)^= %then class &class_switchd &added_class ;;
     model &treatment = &model_switchd   / nologscale   ;
     %if %bquote(&by_statement)^= %then by &by_statement ;; 
     output out =  work._switch_d1 (keep =   &id &period   p1_d ) p = p1_d ;
     run; 

     proc logistic data =  work.sw(where = ( _eligible1 = 1  %if %bquote(&eligible_wts_1)^= %then and &eligible_wts_1 = 1 ; )) descending    &print_option ;
     ods select ModelInfo ResponseProfile ParameterEstimates ;                              
      title3 "Model for P(&treatment = 1 |  &treatment = 1 ,use_censor = &use_censor ) for numerator ";
     %if %bquote(&class_switchn)^= %then class &class_switchn &added_class ;;
     model &treatment = &model_switchn  / nologscale    ;
     %if %bquote(&by_statement)^= %then by &by_statement ;; 
     output out =  work._switch_n1  (keep =  &id &period    p1_n   ) p = p1_n ;
     run; 

     %if %bquote(&cense) ^= %then %do ;
           %if &pool_cense = 1 %then %do;
                 proc logistic data =  work.sw(where = (  &cense > . ))    &print_option ;
                 ods select ModelInfo ResponseProfile ParameterEstimates ;                          
                 title3 "Model for P(&cense = 0 |  X ) for denominator ";
   
                 %if %bquote(&class_censed)^= %then class &class_censed /* &added_class */ ;;
                 model &cense  = &model_censed      / nologscale   ;
                 %if %bquote(&by_statement)^= %then by &by_statement ;; 
                 output out =  work._cense_d0    (keep =  &id &period  pC_d   ) p = pC_d ;
                 run; 
                 
                 proc logistic data =  work.sw(where = (  &cense > . ))  &print_option ;
                 ods select ModelInfo ResponseProfile ParameterEstimates ; 
                 title3 "Model for P(&cense = 0 |  X ) for numerator ";
     
                 %if %bquote(&class_censen)^= %then class &class_censen /* &added_class */ ;;
                 model &cense = &model_censen   / nologscale     ;
                 %if %bquote(&by_statement)^= %then by &by_statement ;; 
                 output out =  work._cense_n0   (keep =  &id &period     pC_n   ) p = pC_n ;
                 run; 
          %end;
          %else %do;

                 proc logistic data =  work.sw(where = ( _eligible0 = 1 AND &cense > . ))    &print_option ;
                 ods select ModelInfo ResponseProfile ParameterEstimates ;                          
                 title3 "Model for P(&cense = 0 |  X, Am1=0 ) for denominator ";
   
                 %if %bquote(&class_censed)^= %then class &class_censed /* &added_class */ ;;
                 model &cense  = &model_censed      / nologscale   ;
                 %if %bquote(&by_statement)^= %then by &by_statement ;; 
                 output out =  work._cense_d0    (keep =  &id &period  pC_d0   ) p = pC_d0 ;
                 run; 

                 proc logistic data =  work.sw(where = (_eligible0 = 1 AND   &cense > . ))  &print_option ;
                 ods select ModelInfo ResponseProfile ParameterEstimates ; 
                 title3 "Model for P(&cense = 0 |  X, Am1=0 ) for numerator ";
     
                 %if %bquote(&class_censen)^= %then class &class_censen /* &added_class */ ;;
                 model &cense = &model_censen   / nologscale     ;
                 %if %bquote(&by_statement)^= %then by &by_statement ;; 
                 output out =  work._cense_n0   (keep =  &id &period     pC_n0   ) p = pC_n0 ;
                 run; 

                 proc logistic data =  work.sw(where = ( _eligible1 = 1 AND &cense > . ))    &print_option ;
                 ods select ModelInfo ResponseProfile ParameterEstimates ;                          
                 title3 "Model for P(&cense = 0 |  X , Am1 = 1) for denominator ";
   
                 %if %bquote(&class_censed)^= %then class &class_censed /* &added_class */ ;;
                 model &cense  = &model_censed      / nologscale   ;
                 %if %bquote(&by_statement)^= %then by &by_statement ;; 
                 output out =  work._cense_d1    (keep =  &id &period  pC_d1   ) p = pC_d1 ;
                 run; 

                 proc logistic data =  work.sw(where = (_eligible1 = 1 AND   &cense > . ))  &print_option ;
                 ods select ModelInfo ResponseProfile ParameterEstimates ; 
                 title3 "Model for P(&cense = 0 |  X, Am1 = 1 ) for numerator ";
     
                 %if %bquote(&class_censen)^= %then class &class_censen /* &added_class */ ;;
                 model &cense = &model_censen   / nologscale     ;
                 %if %bquote(&by_statement)^= %then by &by_statement ;; 
                 output out =  work._cense_n1   (keep =  &id &period     pC_n1   ) p = pC_n1 ;
                 run; 

          %end;

     %end;

                      
     data sw ;
     merge sw work._switch_d0 work._switch_n0 work._switch_d1 work._switch_n1 
                                       %if %bquote(&cense) ^= %then %do;
                                             %if &pool_cense = 1 %then    work._cense_d0 work._cense_n0 ; 
                                             %if &pool_cense = 0 %then    work._cense_d0 work._cense_n0  work._cense_d1 work._cense_n1 ; 
                                    %end;
                               ; /* end of merge */
     by &id &period ;              

 
     /* default setting of weight will be missing for checking below */
     wt = . ;

    /* current observation has _Am1_ = 0 */

    if _Am1_ = 0 then do;
        %if %bquote(&eligible_wts_0)^=  %then %do;
             if &eligible_wts_0 = 1 then do ;
                if      &treatment = 0  and p0_n > . and p0_d > .  then wt = (1.0-p0_n)/(1.0-p0_d);
                else if &treatment = 1  and p0_n > . and p0_d > .  then wt =  p0_n / p0_d ;
             end;
             else if &eligible_wts_0 = 0 then do ;
                wt = 1.0;
             end;
        %end;
        %else %do;

            if      &treatment = 0 and p0_n > . and p0_d > .  then wt = (1.0-p0_n)/(1.0-p0_d);
            else if &treatment = 1 and p0_n > . and p0_d > .  then wt =  p0_n / p0_d ;

        %end;
    end;

   /* current observation has _Am1_ = 1 */
     if _Am1_ = 1 then do ;
         %if %bquote(&eligible_wts_1)^= %then %do;
             if &eligible_wts_1 = 1 then do ;
                 if      &treatment = 0 and p1_n > . and p1_d > .  then wt = (1.0-p1_n)/(1.0-p1_d);
                 else if &treatment = 1 and p1_n > . and p1_d > .  then wt =  p1_n / p1_d ;
             end;
             else if &eligible_wts_1 = 0 then do ;
                  wt = 1 ;
             end;
         %end;
         %else %do;
              if      &treatment = 0 and  p1_n > . and p1_d > .  then wt = (1.0-p1_n)/(1.0-p1_d);
              else if &treatment = 1 and  p1_n > . and p1_d > .  then wt =  p1_n / p1_d ;
         
         %end;

     end;

     %if %bquote(&cense)^= %then %do;
            %if &pool_cense = 0 %then %do;
                if _Am1_ = 0 then do ;
                   pC_n = pC_n0 ;
                   pC_d = pC_d0 ;
                end;
                else if _Am1_ = 1 then do ;
                   pC_n = pC_n1 ;
                   pC_d = pC_d1 ;
                end;
            %end;
            if pC_d = . then pC_d = 1 ;
            if pC_n = . then pC_n = 1 ;                     
            wtC = pC_n/pC_d ;      
     %end;
     %else %do;
        wtC = 1.0 ;
     %end;

     wt = wt * wtC ;

     if wt = . then do ;
        wt = 1.0 ;
        put "problem with weight models "  &id =  &period =  &treatment = _Am1_= p0_d= p0_n= p1_d= p1_n= ;
     end;
     run;
                  
     proc datasets library = work nolist;
     delete  _switch_d0 _switch_d1 _switch_n0 _switch_n1 
             %if %bquote(&cense)^= %then _cense_n0 _cense_d0 ;
           ;
     quit ;

  %end;

    %local ncov i ii varii tmplist word0 word1 ; 

    %let ncov = %numargs(&outcomeCov_var);
    %let tmplist = ;
    %if &ncov > 0 %then %do;
        %let word0 = %scan(&outcomeCov_var,1,%str( ));
        %let tmplist = &word0._new ;
        %do i = 2 %to &ncov;
           %let word0 = %scan(&outcomeCov_var,%eval(&i),%str( ));
           %let tmplist = &tmplist &word0._new ;
        %end;
    %end;
   
    %let array_type = _TEMPORARY_ ;

    data   _switch (keep = &id  _case  _weight_  _assigned_treatment  &tmplist 
                     for_period for_period2 
                     _followup_time _followup_time2 
                      &model_var 
                      )   ;
    set work.sw ;
    by &id &period  ;
    length for_period _followup_time _assigned_treatment 3 for_period2 _followup_time2 5
              
              _case 3
              %do i = 1 %to &ncov ;
                 %let word0 = %scan(&tmplist,%eval(&i),%str( ));
                  &word0 &&length&i 
              %end;
  
     ;
    
    %do i = 1 %to &ncov ;
        %let word0 = %scan(&outcomeCov_var,%eval(&i),%str( ));
        array _&word0 { &first_period:&last_period } &array_type ;       
    %end;

    array elgcount {&first_period:&last_period}   &array_type ;
    array expand   {&first_period:&last_period}   &array_type ;
    array init     {&first_period:&last_period}   &array_type ;
    array wtprod   {&first_period:&last_period}   &array_type ;

    retain _weight0_  ;

    if first.&id then do ;
        do i = &first_period to &last_period ;
            wtprod[i] = 1.0 ;
            elgcount[i] = 0 ;
            expand[i] = 0 ;
             init[i] = . ;
            %do i = 1 %to &ncov ;
               %let word0 = %scan(&outcomeCov_var,%eval(&i),%str( ));
                _&word0 [i] = . ;
            %end;
            
        end;
       _weight0_ = 1.0 ;
    end;

    %if &use_weights = 0 %then wt = 1.0 ;;

    _weight0_ = _weight0_* wt ;

    wtprod[&period] = _weight0_;


    /******************************************************** 
       start basic variables for final model 
       function of treatment                    
     ********************************************************/

     %building_blocks ;
   
    /* end basic variables for final model */

    %if &use_censor = 0 %then switch = 0 ;;

    elgcount[&period] = &eligible ;
    if &eligible = 1 then init[&period] = &assigned_treat ;
    else init[&period] = init[&period - 1] ;

    %do i = 1 %to &ncov ;
        %let word0 = %scan(&outcomeCov_var,%eval(&i),%str( ));
         _&word0 [ &period] = &word0 ;
    %end;

    if &eligible = 1 and &assigned_treat ^= . then expand[&period] = 1 ;

    do for_period = &first_period to &period  ;
        if expand[for_period] = 1 then do ;
            for_period2 = for_period * for_period ;
            _followup_time = &period - for_period + &baseline_offset  ; 
            _followup_time2 = _followup_time ** 2 ;
            _assigned_treatment = init[for_period ];

            /***************************************************************** 
               create variables for final model. First start with dose and dose*dose 
               also include interaction variables for survival curves. Start with two interactions 
             ********************************************************************/
           %create_final_variables ;
            /***************** end of variable definitions ************************/

            if for_period = &period   and (switch = 1 or switch = .)   then do ;
            
              do i = &first_period to &period ; /* at a switch time so no future observations will be
                                       used in expansions prior to current period value */
                 expand[i] = 0 ;
              end;
            end;
            /* weight0 = product of weights upto and including current period */
            /* wtprod = product of weights at a given period */
            %if &lag_p_nosw = 1   %then %do; 
                 /**************************** 
                    when for_period = current period value and we are starting a new trial
                    _weight0_ = wtprod[for_period] and _weight_ = 1.

                    when for_period < current period then this observation will be used for that trial 
                    and the weight should be PROD(m=for_period + 1 to period) wt_m which is
                             Prod(m=1 to period )wt_m / Prod(m=1 to for_period)
                     *******************/     
                   _weight_ = _weight0_ / wtprod[for_period] ;
            %end; 
            %else %do;
                if for_period = &first_period then _weight_ = _weight0_ ;
                else _weight_ = _weight0_ / wtprod[for_period - 1] ;  
            %end;

            %do i = 1 %to &ncov ;
                %let word0 = %scan(&tmplist,%eval(&i),%str( ));
                %let word1 = %scan(&outcomeCov_var,%eval(&i),%str( ));
                &word0 = _&word1 [for_period ];
            %end;
            _case = 0 ;
      
            %if &use_censor = 0 %then %do;
                 /***********************  
                     if we are not censoring then only need to check if outcome at last time. Here
                     we are not censoring due to a switch in treatment 
                 ************************/ 
                
                 if _time_of_event = &period and &outcome = 1 then _case = 1 ;                 
            %end;
            %else %if &use_censor = 1 %then %do;
                /*************************** 
                   when we are censoring due to a switch in treatment, then there are now
                   two cases to consider.
                     1) not the last observation of original followup and there is a switch so
                        the subject is on a different treatment then at the start of the trial.
                        So that all observations in the trial have the same treatment we set _case = .
                        to force proc logistic to drop this observation from model.

                     2) this is the last observation for a subject and switch = 1 so again we set
                        _case = . to drop the observation. if switch = 0 then we will set _case = 1 if
                        an event is 1.
                   ***********************/
                 if switch = 1 then _case = . ;
                 else if switch = 0 and _time_of_event = &period and &outcome = 1 then _case = 1;
            %end;
          /*  newid = newid0  + for_period ;*/
            output _switch ;
          
        end;/* expand[for_period] = 1 */
    end; /* for_period loop */
    run;

  /* expand_new contains all observations used in expanded analysis 
     need to include the baseline variables used in the final model */
 
       proc datasets library = work nolist ;
       modify _switch  ;

       rename 
           _case = &outcome 
           %do i = 1 %to &ncov ;
               %let word0 = %scan(&tmplist,%eval(&i),%str( ));
               %let word1 = %scan(&outcomeCov_var,%eval(&i),%str( ));
               &word0 = &word1 
           %end;
          ; 
       run; 
       quit;

%mend;


%macro building_blocks ;
  /* macro for creating sas arrays for holding previous variables at all previous times used for creating 
     variables used in final model.The time variable is contained in the &period macro variable. This macro 
     should also contain the array definitions.
     The value of array_type is TEMPORARY .

     This macro is used in the cw3 macro in the construction of the final version of the work._switch data set.
  */
        
   
   array _treat_   {&first_period:&last_period}   &array_type ;        
  _treat_[&period]  = &treatment ;

   array _dosesum_ {&first_period:&last_period}   &array_type ;
  _dosesum_[&period] = _cumA_ ;
%mend ;

%macro create_final_variables ;
      /************************************ 
         this macro builds the variables that are used in the final model. This macro should also contain
         any interaction terms that are to be used in the creation of the survival curves. 

        _assigned_treatment holds the assigned treatment variable at the start of the trial 
       ***********************************/
       
        %if &use_censor = 0 %then %do;           
            _dose_ = _cumA_ - _dosesum_[for_period] + _treat_[for_period] ;
            _dose2_ = _dose_ * _dose_ ;
        %end;
        %else %if &use_censor = 1 %then %do;
             &assigned_treat = _assigned_treatment ;
        %end;
  
%mend ;           

%macro for_initial_xbeta;

    %if &use_censor = 0 %then %do;
        _dose_ = 0 ;
        _dose2_ = 0 ;
        _interact1 = 0;
        _interact2 =  0;
    %end;
    %else %if &use_censor = 1 %then %do;
        &assigned_treat = 0 ;
        _interact1 = 0;
       
    %end;
%mend ;

%macro for_untreated ;
     %if &use_censor = 0 %then %do;
         _dose_ = 0 ;
         _dose2_ = 0 ;
         _interact1 = 0 ;
         _interact2 = 0 ;
     %end;
     %else %if &use_censor = 1 %then %do;
         &assigned_treat = 0 ;
         _interact1 = 0 ;
     %end;
%mend ;

%macro for_treated ;
    %if &use_censor = 0 %then %do;
        _dose_ = _followup_time ;
        _dose2_  = _dose_ * _dose_ ;
        _interact1 = _dose_ * _followup_time ;
        _interact2 = _dose2_ * _followup_time;
    %end;
    %else %if &use_censor = 1 %then %do;
         &assigned_treat = 1 ;
         _interact1 = 1 * _followup_time ;
    %end;
%mend ;




/**********************************/


%macro building_blocks0 ;
  /* macro for creating sas arrays for holding previous variables at all previous times used for creating 
     variables used in final model.The time variable is contained in the &period macro variable. This macro 
     should also contain the array definitions.
     The value of array_type is TEMPORARY .

     This macro is used in the cw3 macro in the construction of the final version of the work.sw data set.
  */
        
   
   array _treat_   {&first_period:&last_period}   &array_type ;        
  _treat_[&period]  = &actual_treat ;

   array _dosesum_ {&first_period:&last_period}   &array_type ;
  _dosesum_[&period] = _cumA_ ;
%mend ;

%macro create_final_variables0 ;
      /************************************ 
         this macro builds the variables that are used in the final model. This macro should also contain
         any interaction terms that are to be used in the creation of the survival curves. 
       ***********************************/
      
                
        _dose_ = _cumA_ - _dosesum_[for_period] + _treat_[for_period] ;         
        /****************************
         _dose2_ = _dose_ ** 2 ;
         _dose3_ = _dose_ ** 3 ;
        *******************/

        /** for average cumulative dose **/
        _avg_dose = _dose_ / _followup_time ;

        
%mend ;           

%macro for_initial_xbeta0;

    %if &use_censor = 0 %then %do;
        _avg_dose_ = 0 ;
        
        _interact1 = 0;
      * _interact2 =  0; 
    %end;
    %else %if &use_censor = 1 %then %do;
        &assigned_treat = 0 ;
        _interact1 = 0;
       
    %end;
%mend ;

%macro for_untreated0 ;
     %if &use_censor = 0 %then %do;
       *  _dose_ = 0 ;
       *  _dose3_ = 0 ;
       *  _interact1 = 0 ;
       *  _interact2 = 0 ;
         _avg_dose_ = 0 ;
         _interact1 = 0 ;
     %end;
     %else %if &use_censor = 1 %then %do;
         &assigned_treat = 0 ;
         _interact1 = 0 ;
     %end;
%mend ;

%macro for_treated0 ;
    %if &use_censor = 0 %then %do;
       * _dose_ = _followup_time ;
       * _dose3_  = _dose_ * _dose_ * _dose_;
       * _interact1 = _dose_ * _followup_time ;
       * _interact2 = _dose3_ * _followup_time;
         _avg_dose_ = 1 ;
         _interact1 = 1 * _followup_time ;
    %end;
    %else %if &use_censor = 1 %then %do;
         &assigned_treat = 1 ;
         _interact1 = 1 * _followup_time ;
    %end;
%mend ;









