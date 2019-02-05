* Feel free to report errors/suggestions to ehsan@stat.ubc.ca

/***************************************************************
PROGRAM 13.1
Estimating the mean outcome within levels of treatment 
and confounders: Data from NHEFS
***************************************************************/

* some preprocessing of the data
clear
insheet using "C:\Ehsan\nhefs.csv"

gen byte cens = (wt82 == .)
gen byte older = (age > 50 & age < .)
label var older "50+ years"
recode school (0/8 = 1 "1. 8th grade or less") (9/11 = 2 "2. HS dropout")(12 = 3 "3. HS") (13/15 = 4 "4. College dropout") (16/max = 5 "5. College or more") (missing = 6 "Unknown"), gen(education)

* provisionally ignore subjects with missing values
* Analysis restricted to N=1566
foreach var of varlist qsmk sex race age school smokeintensity smokeyrs exercise active wt71 wt82 {
drop if `var'==.
}

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

*Table
foreach var of varlist years male white university kg cigs noexer inactive {
tabdisp qsmk, cell(`var') format(%3.1f)
}

* Estimates 
glm wt82_71 qsmk sex race c.age##c.age ib(last).education c.smokeintensity##c.smokeintensity c.smokeyrs##c.smokeyrs ib(last).exercise ib(last).active c.wt71##c.wt71 qsmk##c.smokeintensity
predict meanY
summarize meanY

list meanY if seqn == 24770

* for comparison with observed outcome
summarize wt82_71


/***************************************************************
PROGRAM 13.2
Standardizing the mean outcome to the baseline confounders
Data from Table 2.2
***************************************************************/
ssc inst tomata
clear
clear mata
mata
obsID = ("Rheia", "Kronos", "Demeter", "Hades", "Hestia", "Poseidon", "Hera", "Zeus", "Artemis", "Apollo", "Leto", "Ares", "Athena", "Hephaestus", "Aphrodite", "Cyclope", "Persephone", "Hermes", "Hebe", "Dionysus")'
obsL = (0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)'
obsA = (0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1)'
obsY = (0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0)'
N = length(obsID)
block = J(1, N, -1)' \ J(1, N, 0)' \ J(1, N, 1)'
L = obsL \ obsL \ obsL
A = obsA \ J(1, N, 0)' \ J(1, N, 1)'
Y = obsY \ J(1, N, .)' \ J(1, N, .)'
id = obsID \ obsID \ obsID
end

getmata Y L A id block
glm Y A##L
predict meanY
summarize meanY

mean meanY if block == -1
mean meanY if block == 0
mean meanY if block == 1

drop meanY

* creation of bootstrap program
capture program drop gboottest
program define gboottest, rclass
preserve
bsample
glm Y A##L
predict meanY, mu 
summarize meanY if block == -1, meanonly
return scalar mtxObs = r(mean)
summarize meanY if block == 1, meanonly
return scalar mtx1 = r(mean)
summarize meanY if block == 0, meanonly
return scalar mtx0 = r(mean)
return scalar mtxD = return(mtx1) - return(mtx0)
drop meanY
restore
end

* number of bootstrap samples to run is set to 50
set seed 1232 
bootstrap r(mtxObs) r(mtx1) r(mtx0) r(mtxD), reps(50): gboottest
estat boot, all


/***************************************************************
 PROGRAM 13.3
Standardizing the mean outcome to the baseline confounders:
Data from NHEFS
***************************************************************/

* Re-run the PROGRAM 13.1 (see above) and then continue as follows:
* create a dataset with 3 copies of each subject
* 1st copy: equal to original one
* 2nd copy: treatment set to 0, outcome to missing
* 3rd copy: treatment set to 1, outcome to missing

putmata *

mata
N = length(wt82_71)
block = J(1, N, -1)' \ J(1, N, 0)' \ J(1, N, 1)'

wt82_71 = wt82_71 \ J(1, N, .)' \ J(1, N, .)'
qsmk = qsmk \ J(1, N, 0)' \ J(1, N, 1)'
sex = sex  \ sex \ sex
race = race \ race \ race
age = age \ age \ age
education = education \ education \ education
smokeintensity = smokeintensity \ smokeintensity \ smokeintensity
smokeyrs = smokeyrs \ smokeyrs \ smokeyrs
exercise = exercise \ exercise \ exercise
active = active \ active \ active
wt71 = wt71 \ wt71 \ wt71
end

clear
getmata block wt82_71 qsmk sex race age education smokeintensity smokeyrs exercise active wt71

* linear model to estimate mean outcome conditional on treatment and confounders;
* parameters are estimated using original observations only (block= -1) ;
* parameter estimates are used to predict mean outcome for observations with 
* treatment set to 0 (block=0) and to 1 (block=1);

glm wt82_71 qsmk sex race c.age##c.age ib(last).education c.smokeintensity##c.smokeintensity c.smokeyrs##c.smokeyrs ib(last).exercise ib(last).active c.wt71##c.wt71
predict meanY
summarize meanY
summarize wt82_71

* estimate mean outcome in each of the groups block=0, and block=1;
* this mean outcome is a weighted average of the mean outcomes in each combination 
* of values of treatment and confounders, that is, the standardized outcome;

mean meanY if block == -1
mean meanY if block == 0
mean meanY if block == 1


/***************************************************************
PROGRAM 13.4
Computing the 95% confidence interval of the standardized means 
and their difference: Data from NHEFS
***************************************************************/

* Re-run the PROGRAM 13.3 (see above) and then continue as follows:
drop  meanY

* creation of bootstrap program
capture program drop gboot
program define gboot, rclass
preserve
bsample
glm wt82_71 qsmk sex race c.age##c.age ib(last).education c.smokeintensity##c.smokeintensity c.smokeyrs##c.smokeyrs ib(last).exercise ib(last).active c.wt71##c.wt71
predict meanY, mu 
summarize meanY if block == -1, meanonly
return scalar mtxObs = r(mean)
summarize meanY if block == 1, meanonly
return scalar mtx1 = r(mean)
summarize meanY if block == 0, meanonly
return scalar mtx0 = r(mean)
return scalar mtxD = return(mtx1) - return(mtx0)
drop meanY
restore
end

* number of bootstrap samples to run is set to 1000 
* to test the program use 50 or so
set seed 1232
bootstrap r(mtxObs) r(mtx1) r(mtx0) r(mtxD), reps(1000): gboot
estat boot, all



