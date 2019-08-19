/* https://www.hsph.harvard.edu/causal/software/ */
/* https://cdn1.sph.harvard.edu/wp-content/uploads/sites/148/2012/10/samplemsmcox.txt */
/*
 	An earlier version of this sample SAS program appeared in the Appendix of

	Hernán MA, Brumback B, Robins JM. "Marginal structural models 
		to estimate the causal effect of zidovudine on the survival 
		of HIV-positive men" 
		Epidemiology 2000;11(5):561-570

	Please refer to that paper for an explanation of the program
*/


/* Model 1 */
proc logistic data=MAIN;
	where MONTH<=ZDV_M or ZDV_M=.;
        model A = AGE_0 YEAR_01 YEAR_02 YEAR_03 CD4_01 CD4_02 CD8_01 CD8_02 
			WBC_01 WBC_02 RBC_01 RBC_02 PLAT_01 PLAT_02 SYMPT_0 
			MONTH MONTH1-MONTH3;
        output out=model1 p=pzdv_0;
run;

/* Model 2 */
proc logistic data=MAIN;
	where MONTH<=ZDV_m or ZDV_M=.;
        model A = AGE_0 YEAR_01 YEAR_02 YEAR_03 CD4_01 CD4_02 CD8_01 CD8_02 
			WBC_01 WBC_02 RBC_01 RBC_02 PLATE_01 PLATE_02 SYMPT_0 
			CD4_1 CD4_2 CD8_1 CD8_2 WBC_1 WBC_2 RBC_1 RBC_2 PLAT_1 PLAT_2
			SYMPT AIDS MONTH MONTH1-MONTH3;
        output out=model2 p=pzdv_w;
run;

/* Model 3 */
proc logistic data=MAIN;
        model C = A AGE_0 YEAR_01 YEAR_02 YEAR_03 CD4_01 CD4_02 CD8_01 CD8_02 
			WBC_01 WBC_02 RBC_01 RBC_02 PLATE_01 PLATE_02 SYMPT_0 
			MONTH MONTH1-MONTH3;
 output out=model3 p=punc_0;
run;

/* Model 4 */
proc logistic data=MAIN;
        model C = A AGE_0 YEAR_01 YEAR_02 YEAR_03 CD4_01 CD4_02 CD8_01 CD8_02 
			WBC_01 WBC_02 RBC_01 RBC_02 PLATE_01 PLATE_02 SYMPT_0 
			CD4_1 CD4_2 CD8_1 CD8_2 WBC_1 WBC_2 RBC_1 RBC_2 PLAT_1 PLAT_2
			SYMPT AIDS MONTH MONTH1-MONTH3;
 output out=model4 p=punc_w;
run;

data main_w;
        merge model1 model2 model3 model4;
        by ID MONTH;
/* variables ending with _0 refer to the numerator of the weights
   Variables ending with _w refer to the denominator of the weights */

/* reset the variables for a new patient */
        if first.id then do; 
	k1_0=1; k2_0=1; k1_w=1; k2_w=1; 
	end;
        retain k1_0 k2_0 k1_w k2_w;

/* Inverse probability of censoring weights */
        k2_0=k2_0*punc_0;
        k2_w=k2_w*punc_w;

/* Inverse probability of treatment weights */
	 /* patients not on zidovudine */
        if zdv_m>month or zdv_m=. then do;
                        	k1_0=k1_0*pzdv_0;
				k1_w=k1_w*pzdv_w;
                        end;
	 /* patients that start zidovudine in the current month */
                else if zdv_m=month then do;
                        	k1_0=k1_0*(1-pzdv_0); 
				k1_w=k1_w*(1-pzdv_w);
                        end;
	 /* patients that have already started zidovudine */
                else do; 
				k1_0=k1_0;
				k1_w=k1_w; 
			    end;

/* Stabilized and non stabilized weights */
        stabw=(k1_0*k2_0)/(k1_w*k2_w);
        nstabw=1/(k1_w*k2_w);
run;

proc genmod data=main_w;
        class id;
        model D= A AGE_0 YEAR_01 YEAR_02 YEAR_03 CD4_01 CD4_02 CD8_01 CD8_02 
			WBC_01 WBC_02 RBC_01 RBC_02 PLAT_01 PLAT_02 SYMPT_0 
			MONTH MONTH1-MONTH3/ link=logit dist=bin;
        weight stabw;
        repeated subject=ID/ type=ind;
run;
