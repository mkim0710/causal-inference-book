** Programmed by Min-Hyung Kim
** Please find the most updated program at the following URL:
** https://github.com/mkim0710/causal-inference-book/blob/master/stata/standardization_parametric_regress.do


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


     +-------------------------------------------------------------------------------------------------------------+
     | qsmk   sex   race   exercise   age_gt50     n   _fillin   E_outcome   total   n_by_c~r   E_outc~x      diff |
     |-------------------------------------------------------------------------------------------------------------|
  1. |    0     0      0          0          0   101         0    3.441655    1566        130   1.871614   2.97923 |
  2. |    0     0      0          1          0   142         0    3.211052    1566        201   1.871614   2.97923 |
  3. |    0     0      0          2          0    93         0    2.994104    1566        128   1.871614   2.97923 |
  4. |    0     0      1          0          0    10         0    3.364141    1566         11   1.871614   2.97923 |
  5. |    0     0      1          1          0    15         0    3.133537    1566         18   1.871614   2.97923 |
     |-------------------------------------------------------------------------------------------------------------|
  6. |    0     0      1          2          0    22         0    2.916589    1566         27   1.871614   2.97923 |
  7. |    0     1      0          0          0    77         0    3.136499    1566         89   1.871614   2.97923 |
  8. |    0     1      0          1          0   173         0    2.905895    1566        221   1.871614   2.97923 |
  9. |    0     1      0          2          0   142         0    2.688948    1566        182   1.871614   2.97923 |
 10. |    0     1      1          0          0     3         0    3.058985    1566          4   1.871614   2.97923 |
     |-------------------------------------------------------------------------------------------------------------|
 11. |    0     1      1          1          0    26         0    2.828381    1566         34   1.871614   2.97923 |
 12. |    0     1      1          2          0    47         0    2.611433    1566         53   1.871614   2.97923 |
 13. |    0     0      0          0          1    34         0    -.349313    1566         43   1.871614   2.97923 |
 14. |    0     0      0          1          1    52         0   -.5799167    1566         87   1.871614   2.97923 |
 15. |    0     0      0          2          1    47         0   -.7968646    1566         83   1.871614   2.97923 |
     |-------------------------------------------------------------------------------------------------------------|
 16. |    0     0      1          0          1     5         0   -.4268273    1566          6   1.871614   2.97923 |
 17. |    0     0      1          1          1    10         0    -.657431    1566         13   1.871614   2.97923 |
 18. |    0     0      1          2          1    11         0   -.8743789    1566         15   1.871614   2.97923 |
 19. |    0     1      0          0          1     7         0    -.654469    1566         17   1.871614   2.97923 |
 20. |    0     1      0          1          1    63         0   -.8850727    1566         83   1.871614   2.97923 |
     |-------------------------------------------------------------------------------------------------------------|
 21. |    0     1      0          2          1    62         0   -1.102021    1566         96   1.871614   2.97923 |
 22. |    0     1      1          0          1     .         1   -.7319834    1566          0   1.871614   2.97923 |
 23. |    0     1      1          1          1     4         0   -.9625871    1566          4   1.871614   2.97923 |
 24. |    0     1      1          2          1    17         0   -1.179535    1566         21   1.871614   2.97923 |
 25. |    1     0      0          0          0    29         0    6.420885    1566        130   4.850844   2.97923 |
     |-------------------------------------------------------------------------------------------------------------|
 26. |    1     0      0          1          0    59         0    6.190281    1566        201   4.850844   2.97923 |
 27. |    1     0      0          2          0    35         0    5.973333    1566        128   4.850844   2.97923 |
 28. |    1     0      1          0          0     1         0    6.343371    1566         11   4.850844   2.97923 |
 29. |    1     0      1          1          0     3         0    6.112767    1566         18   4.850844   2.97923 |
 30. |    1     0      1          2          0     5         0    5.895819    1566         27   4.850844   2.97923 |
     |-------------------------------------------------------------------------------------------------------------|
 31. |    1     1      0          0          0    12         0    6.115729    1566         89   4.850844   2.97923 |
 32. |    1     1      0          1          0    48         0    5.885125    1566        221   4.850844   2.97923 |
 33. |    1     1      0          2          0    40         0    5.668178    1566        182   4.850844   2.97923 |
 34. |    1     1      1          0          0     1         0    6.038215    1566          4   4.850844   2.97923 |
 35. |    1     1      1          1          0     8         0    5.807611    1566         34   4.850844   2.97923 |
     |-------------------------------------------------------------------------------------------------------------|
 36. |    1     1      1          2          0     6         0    5.590663    1566         53   4.850844   2.97923 |
 37. |    1     0      0          0          1     9         0    2.629917    1566         43   4.850844   2.97923 |
 38. |    1     0      0          1          1    35         0    2.399313    1566         87   4.850844   2.97923 |
 39. |    1     0      0          2          1    36         0    2.182365    1566         83   4.850844   2.97923 |
 40. |    1     0      1          0          1     1         0    2.552402    1566          6   4.850844   2.97923 |
     |-------------------------------------------------------------------------------------------------------------|
 41. |    1     0      1          1          1     3         0    2.321799    1566         13   4.850844   2.97923 |
 42. |    1     0      1          2          1     4         0    2.104851    1566         15   4.850844   2.97923 |
 43. |    1     1      0          0          1    10         0    2.324761    1566         17   4.850844   2.97923 |
 44. |    1     1      0          1          1    20         0    2.094157    1566         83   4.850844   2.97923 |
 45. |    1     1      0          2          1    34         0    1.877209    1566         96   4.850844   2.97923 |
     |-------------------------------------------------------------------------------------------------------------|
 46. |    1     1      1          0          1     .         1    2.247246    1566          0   4.850844   2.97923 |
 47. |    1     1      1          1          1     .         1    2.016643    1566          4   4.850844   2.97923 |
 48. |    1     1      1          2          1     4         0    1.799695    1566         21   4.850844   2.97923 |
     +-------------------------------------------------------------------------------------------------------------+

*/
