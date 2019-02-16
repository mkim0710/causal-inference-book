clear
use "D:\OneDrive - Harvard University\[[[[School]]]]\[[[Harvard MPH]]]\[[EPI524]]\nhefs.dta" 
gen age_gt50 = (age > 50 & age <.)
gen age_ge50 = (age >= 50 & age <.)
drop if wt82_71 == .
global varname4id seqn
global varname4outcome wt82_71
global varname4tx qsmk
global varnames4confounder age_gt50
restore, preserve
preserve
collapse (count) n = $varname4id (mean) E_outcome = $varname4outcome, by($varname4tx $varnames4confounder)
egen total = sum(n)
gen p = n/total
egen n_by_confounder = sum(n), by($varnames4confounder)
egen p_by_confounder = sum(p), by($varnames4confounder)
list
drop p p_by_confounder
egen E_outcome_by_tx = sum(E_outcome * n_by_confounder / total), by($varname4tx)
egen tmp_min = min(E_outcome_by_tx)
egen tmp_max = max(E_outcome_by_tx)
gen diff = tmp_max - tmp_min
drop tmp*
list
restore, preserve
/*

     +----------------------------------------------------------------------------+
     | qsmk   older     n      E_outcome   total   n_by_c~r   E_outc~x       diff |
     |----------------------------------------------------------------------------|
  1. |    0       0   851    2.991448404    1566       1098   1.869721   3.004455 |
  2. |    0       1   312   -.7620255091    1566        468   1.869721   3.004455 |
  3. |    1       0   247    6.056408169    1566       1098   4.874175   3.004455 |
  4. |    1       1   156    2.100474456    1566        468   4.874175   3.004455 |
     +----------------------------------------------------------------------------+
*/
