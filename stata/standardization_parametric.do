** Programmed by Min-Hyung Kim
** Please find the most updated program at the following URL:
** https://github.com/mkim0710/causal-inference-book/blob/master/stata/standardization_parametric.do


** The dataset nhefs.dta is from the following reference:  
** Hernán MA, Robins JM (2019). Causal Inference. Boca Raton: Chapman & Hall/CRC, forthcoming. 
** https://www.hsph.harvard.edu/miguel-hernan/causal-inference-book/


clear
webuse set "https://github.com/mkim0710/causal-inference-book/raw/master/data"
webuse "nhefs.dta"
gen age_gt50 = (age > 50 & age <.)
drop if wt82_71 == .
global varname4id seqn
global varname4outcome wt82_71
global varname4tx qsmk
global varname4confounder age_gt50
global i_varname4confounder
foreach varname of global varname4confounder {
    global i_varname4confounder $i_varname4confounder i.`varname'
}
regress $varname4outcome $varname4tx $i_varname4confounder
preserve
collapse (count) n = $varname4id, by($varname4tx $varname4confounder)
fillin $varname4tx $varname4confounder
predict E_outcome, xb
egen total = sum(n)
egen n_by_confounder = sum(n), by($varname4confounder)
egen E_outcome_by_tx = sum(E_outcome * n_by_confounder / total), by($varname4tx)
egen tmp_min = min(E_outcome_by_tx)
egen tmp_max = max(E_outcome_by_tx)
gen diff = tmp_max - tmp_min
drop tmp*
list
restore, preserve

/*
     +---------------------------------------------------------------------------+
     | qsmk   age_gt50     n   E_outcome   total   n_by_c~r   E_outc~x      diff |
     |---------------------------------------------------------------------------|
  1. |    0          0   851    3.007481    1566       1098   1.867893   2.99369 |
  2. |    0          1   312   -.8057554    1566        468   1.867893   2.99369 |
  3. |    1          0   247    6.001171    1566       1098   4.861583   2.99369 |
  4. |    1          1   156    2.187934    1566        468   4.861583   2.99369 |
     +---------------------------------------------------------------------------+
*/

global varname4id seqn
global varname4outcome wt82_71
global varname4tx qsmk
global varname4confounder age_gt50 sex race exercise
global i_varname4confounder
foreach varname of global varname4confounder {
    global i_varname4confounder $i_varname4confounder i.`varname'
}
regress $varname4outcome $varname4tx $i_varname4confounder
preserve
collapse (count) n = $varname4id, by($varname4tx $varname4confounder)
fillin $varname4tx $varname4confounder
predict E_outcome, xb
egen total = sum(n)
egen n_by_confounder = sum(n), by($varname4confounder)
egen E_outcome_by_tx = sum(E_outcome * n_by_confounder / total), by($varname4tx)
egen tmp_min = min(E_outcome_by_tx)
egen tmp_max = max(E_outcome_by_tx)
gen diff = tmp_max - tmp_min
drop tmp*
list
restore, preserve

/*

     +----------------------------------------------------------------------------------------------------+
     | qsmk   sex   race   exercise   age_gt50     n   E_outcome   total   n_by_c~r   E_outc~x       diff |
     |----------------------------------------------------------------------------------------------------|
  1. |    0     0      0          0          0   101    3.437669    1566        130   1.871671   2.973847 |
  2. |    0     0      0          1          0   142    3.214942    1566        201   1.871671   2.973847 |
  3. |    0     0      0          2          0    93    2.992215    1566        128   1.871671   2.973847 |
  4. |    0     0      1          0          0    10    3.360514    1566         11   1.871671   2.973847 |
  5. |    0     0      1          1          0    15    3.137787    1566         18   1.871671   2.973847 |
     |----------------------------------------------------------------------------------------------------|
  6. |    0     0      1          2          0    22     2.91506    1566         27   1.871671   2.973847 |
  7. |    0     1      0          0          0    77    3.132217    1566         89   1.871671   2.973847 |
  8. |    0     1      0          1          0   173     2.90949    1566        221   1.871671   2.973847 |
  9. |    0     1      0          2          0   142    2.686763    1566        182   1.871671   2.973847 |
 10. |    0     1      1          0          0     3    3.055062    1566          4   1.871671   2.973847 |
     |----------------------------------------------------------------------------------------------------|
 11. |    0     1      1          1          0    26    2.832335    1566         34   1.871671   2.973847 |
 12. |    0     1      1          2          0    47    2.609608    1566         53   1.871671   2.973847 |
 13. |    0     0      0          0          1    34   -.3532559    1566         43   1.871671   2.973847 |
 14. |    0     0      0          1          1    52   -.5759827    1566         87   1.871671   2.973847 |
 15. |    0     0      0          2          1    47   -.7987096    1566         83   1.871671   2.973847 |
     |----------------------------------------------------------------------------------------------------|
 16. |    0     0      1          0          1     5   -.4304104    1566          6   1.871671   2.973847 |
 17. |    0     0      1          1          1    10   -.6531372    1566         13   1.871671   2.973847 |
 18. |    0     0      1          2          1    11   -.8758641    1566         15   1.871671   2.973847 |
 19. |    0     1      0          0          1     7   -.6587077    1566         17   1.871671   2.973847 |
 20. |    0     1      0          1          1    63   -.8814346    1566         83   1.871671   2.973847 |
     |----------------------------------------------------------------------------------------------------|
 21. |    0     1      0          2          1    62   -1.104161    1566         96   1.871671   2.973847 |
 22. |    0     1      1          1          1     4   -.9585891    1566          4   1.871671   2.973847 |
 23. |    0     1      1          2          1    17   -1.181316    1566         21   1.871671   2.973847 |
 24. |    1     0      0          0          0    29    6.416677    1566        130   4.845519   2.973847 |
 25. |    1     0      0          1          0    59     6.19395    1566        201   4.845519   2.973847 |
     |----------------------------------------------------------------------------------------------------|
 26. |    1     0      0          2          0    35    5.971223    1566        128   4.845519   2.973847 |
 27. |    1     0      1          0          0     1    6.339522    1566         11   4.845519   2.973847 |
 28. |    1     0      1          1          0     3    6.116795    1566         18   4.845519   2.973847 |
 29. |    1     0      1          2          0     5    5.894068    1566         27   4.845519   2.973847 |
 30. |    1     1      0          0          0    12    6.111225    1566         89   4.845519   2.973847 |
     |----------------------------------------------------------------------------------------------------|
 31. |    1     1      0          1          0    48    5.888498    1566        221   4.845519   2.973847 |
 32. |    1     1      0          2          0    40    5.665771    1566        182   4.845519   2.973847 |
 33. |    1     1      1          0          0     1     6.03407    1566          4   4.845519   2.973847 |
 34. |    1     1      1          1          0     8    5.811343    1566         34   4.845519   2.973847 |
 35. |    1     1      1          2          0     6    5.588616    1566         53   4.845519   2.973847 |
     |----------------------------------------------------------------------------------------------------|
 36. |    1     0      0          0          1     9    2.625752    1566         43   4.845519   2.973847 |
 37. |    1     0      0          1          1    35    2.403025    1566         87   4.845519   2.973847 |
 38. |    1     0      0          2          1    36    2.180298    1566         83   4.845519   2.973847 |
 39. |    1     0      1          0          1     1    2.548598    1566          6   4.845519   2.973847 |
 40. |    1     0      1          1          1     3    2.325871    1566         13   4.845519   2.973847 |
     |----------------------------------------------------------------------------------------------------|
 41. |    1     0      1          2          1     4    2.103144    1566         15   4.845519   2.973847 |
 42. |    1     1      0          0          1    10      2.3203    1566         17   4.845519   2.973847 |
 43. |    1     1      0          1          1    20    2.097573    1566         83   4.845519   2.973847 |
 44. |    1     1      0          2          1    34    1.874847    1566         96   4.845519   2.973847 |
 45. |    1     1      1          2          1     4    1.797692    1566         21   4.845519   2.973847 |
     +----------------------------------------------------------------------------------------------------+
*/

