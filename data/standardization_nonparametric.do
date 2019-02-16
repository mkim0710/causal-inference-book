clear
webuse set "https://github.com/mkim0710/causal-inference-book/raw/master/data"
webuse "nhefs.dta"
gen age_gt50 = (age > 50 & age <.)
drop if wt82_71 == .
global varname4id seqn
global varname4outcome wt82_71
global varname4tx qsmk
global varnames4confounder age_gt50
preserve
collapse (count) n = $varname4id (mean) E_outcome = $varname4outcome, by($varname4tx $varnames4confounder)
egen total = sum(n)
egen n_by_confounder = sum(n), by($varnames4confounder)
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


global varname4id seqn
global varname4outcome wt82_71
global varname4tx qsmk
global varnames4confounder age_gt50 sex race
preserve
collapse (count) n = $varname4id (mean) E_outcome = $varname4outcome, by($varname4tx $varnames4confounder)
egen total = sum(n)
egen n_by_confounder = sum(n), by($varnames4confounder)
egen E_outcome_by_tx = sum(E_outcome * n_by_confounder / total), by($varname4tx)
egen tmp_min = min(E_outcome_by_tx)
egen tmp_max = max(E_outcome_by_tx)
gen diff = tmp_max - tmp_min
drop tmp*
list
restore, preserve
/*

     +--------------------------------------------------------------------------------------------+
     | qsmk   sex   race   age_gt50     n      E_outcome   total   n_by_c~r   E_outc~x       diff |
     |--------------------------------------------------------------------------------------------|
  1. |    0     0      0          0   336    3.068112919    1566        459   1.885048   2.989476 |
  2. |    0     0      1          0    47    2.770156591    1566         56   1.885048   2.989476 |
  3. |    0     1      0          0   392    2.883348802    1566        492   1.885048   2.989476 |
  4. |    0     1      1          0    76    3.346928431    1566         91   1.885048   2.989476 |
  5. |    0     0      0          1   133   -.2280019062    1566        213   1.885048   2.989476 |
     |--------------------------------------------------------------------------------------------|
  6. |    0     0      1          1    26    -1.77060433    1566         34   1.885048   2.989476 |
  7. |    0     1      0          1   132   -.8669447827    1566        196   1.885048   2.989476 |
  8. |    0     1      1          1    21   -2.235965782    1566         25   1.885048   2.989476 |
  9. |    1     0      0          0   123     6.15771731    1566        459   4.874525   2.989476 |
 10. |    1     0      1          0     9    8.292129953    1566         56   4.874525   2.989476 |
     |--------------------------------------------------------------------------------------------|
 11. |    1     1      0          0   100    5.940416326    1566        492   4.874525   2.989476 |
 12. |    1     1      1          0    15    4.657519101    1566         91   4.874525   2.989476 |
 13. |    1     0      0          1    80    2.624831565    1566        213   4.874525   2.989476 |
 14. |    1     0      1          1     8    2.636879605    1566         34   4.874525   2.989476 |
 15. |    1     1      0          1    64    1.255936479    1566        196   4.874525   2.989476 |
     |--------------------------------------------------------------------------------------------|
 16. |    1     1      1          1     4    4.053129605    1566         25   4.874525   2.989476 |
     +--------------------------------------------------------------------------------------------+
*/
