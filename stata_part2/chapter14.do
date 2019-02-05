* Feel free to report errors/suggestions to ehsan@stat.ubc.ca

/***************************************************************
# PROGRAM 14.1
# Preprocessing, ranks of extreme observations, IP weights for censoring
# Data from NHEFS
***************************************************************/

* some preprocessing of the data
clear
insheet using "C:\Ehsan\nhefs.csv"

gen byte cens = (wt82 == .)
gen byte older = (age > 50 & age < .)
label var older "50+ years"
recode school (0/8 = 1 "1. 8th grade or less") (9/11 = 2 "2. HS dropout")(12 = 3 "3. HS") (13/15 = 4 "4. College dropout") (16/max = 5 "5. College or more") (missing = 6 "Unknown"), gen(education)

* provisionally ignore subjects with missing values in school variable
drop if school ==.

* mean weight change in those with and without smoking cessation
label define qsmk 0 "No smoking cessation" 1 "Smoking cessation"
label values qsmk qsmk
label var wt82_71 "weight change in kg"
bysort qsmk: sum wt82_71 if cens == 0
regress wt82_71 qsmk if cens == 0, noheader
by qsmk, sort: egen years = mean(age) if age < . & cens == 0
label var years "Age, years"
by qsmk, sort: egen cigs = mean(smokeintensity) if smokeintensity < . & cens== 0
label var cigs "Cigarettes/day"
by qsmk, sort: egen kg = mean(wt71) if wt71 < . & cens == 0
label var kg "Weight, kg"
by qsmk, sort: egen white = mean(100 * (race==0)) if race < . & cens == 0
label var white "White, %"
by qsmk, sort: egen male = mean(100 * (sex==0)) if sex < . & cens == 0
label var male "Men, %"
by qsmk, sort: egen inactive = mean(100 * (active==2)) if active < . & cens ==0
label var inactive "Inactive daily life"
by qsmk, sort: egen noexer = mean(100 * (exercise == 2)) if exercise < . & cens == 0
label var noexer "Little/no exercise"
by qsmk, sort: egen university = mean(100 * (education == 5)) if education <. & cens == 0
label var university "University education"

* estimation of denominator of ip weights for C
glm cens qsmk sex race c.age##c.age ib(last).education c.smokeintensity##c.smokeintensity c.smokeyrs##c.smokeyrs ib(last).exercise ib(last).active c.wt71##c.wt71
predict pr_cens
gen w_cens = 1/(1-pr_cens)
summarize w_cens

* Analysis restricted to N=1566
drop if wt82 == .
summarize wt82_71


/***************************************************************
#  PROGRAM 14.2
# G-estimation of a 1-parameter structural nested mean model
# Brute force search
# Data from NHEFS
***************************************************************/

gen psi = 3.446
gen Hpsi = wt82_71 - psi * qsmk 

xi: xtgee qsmk sex race c.age##c.age ib(last).education c.smokeintensity##c.smokeintensity c.smokeyrs##c.smokeyrs ib(last).exercise ib(last).active c.wt71##c.wt71 Hpsi [pw = w_cens], link(logit) i(seqn) corr(ind) family(binomial)
di _b[Hpsi]

drop psi Hpsi

/***************************************************************
# G-estimation: Checking multiple possible values of psi
***************************************************************/

local seq_start = 2
local seq_end = 5
local seq_by = 0.1
local seq_len = (`seq_end'-`seq_start')/`seq_by' + 1
matrix results = J(`seq_len',3,0)
matrix list results
gen psi = .
gen Hpsi = .
local j = 0
forvalues i =  `seq_start'(`seq_by')`seq_end' {
local j = `j'+1
replace psi = `i'
replace Hpsi = wt82_71 - psi * qsmk 
quietly xtgee qsmk sex race c.age##c.age ib(last).education c.smokeintensity##c.smokeintensity c.smokeyrs##c.smokeyrs ib(last).exercise ib(last).active c.wt71##c.wt71 Hpsi [pw = w_cens], link(logit) i(seqn) corr(ind) family(binomial)
local b = _b[Hpsi]
di "coeff `b' is generated from psi `i'"
matrix results[`j',1]= `i'
matrix results[`j',2]= `b'
matrix results[`j',3]= abs(`b')
}
matrix colnames results = "psi" "B(Hpsi)" "AbsB(Hpsi)"
mat li results 

mata
res = st_matrix("results")
for(i=1; i<= rows(res); i++) { 
if (res[i,3] == colmin(res[,3])) res[i,1]	
}
end
* Setting seq_by = 0.01 will yield the result 3.46

/***************************************************************
#  PROGRAM 14.3
# G-estimation for 2-parameter structural nested mean model
# Closed form estimator
# Data from NHEFS
***************************************************************/

/***************************************************************
# G-estimation: Closed form estimator linear mean models  
***************************************************************/
  
glm qsmk sex race c.age##c.age ib(last).education c.smokeintensity##c.smokeintensity c.smokeyrs##c.smokeyrs ib(last).exercise ib(last).active c.wt71##c.wt71[pw = w_cens], link(logit) family(binomial)
predict pr_qsmk
summarize pr_qsmk

ssc inst tomata
putmata *
mata: diff = qsmk - pr_qsmk
mata: num = w_cens :* wt82_71 :* diff
mata: den = sum(w_cens :* qsmk :* diff)
mata: sum( num/den )

/***************************************************************
# G-estimation: Closed form estimator for 2-parameter model
***************************************************************/
mata
diff = qsmk - pr_qsmk
diff2 = w_cens :* diff

lhs = J(2,2, 0)
lhs[1,1] = sum( qsmk :* diff2)
lhs[1,2] = sum( qsmk :* smokeintensity  :* diff2 )
lhs[2,1] = sum( qsmk :* smokeintensity :* diff2)
lhs[2,2] = sum( qsmk :* smokeintensity :* smokeintensity :* diff2 )
                                                                
rhs = J(2,1,0)
rhs[1] = sum(wt82_71 :* diff2 )
rhs[2] = sum(wt82_71 :* smokeintensity :* diff2 )

psi = (lusolve(lhs, rhs))'
psi
psi = (invsym(lhs'lhs)*lhs'rhs)'
psi
end



