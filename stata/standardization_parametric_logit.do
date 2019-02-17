** Programmed by Min-Hyung Kim
** Please find the most updated program at the following URL:
** https://github.com/mkim0710/causal-inference-book/blob/master/stata/standardization_parametric_logit.do


** The dataset nhefs.dta is from the following reference:  
** Hernán MA, Robins JM (2019). Causal Inference. Boca Raton: Chapman & Hall/CRC, forthcoming. 
** https://www.hsph.harvard.edu/miguel-hernan/causal-inference-book/


clear
webuse set "https://github.com/mkim0710/causal-inference-book/raw/master/data"
webuse "nhefs.dta"
gen age_gt50 = (age > 50 & age <.)
global varname4id seqn
global varname4outcome death
global varname4tx qsmk
global varname4confounder age_gt50
global i_varname4confounder
foreach varname of global varname4confounder {
    global i_varname4confounder $i_varname4confounder i.`varname'
}
logit $varname4outcome $varname4tx $i_varname4confounder
preserve
collapse (count) n = $varname4id, by($varname4tx $varname4confounder)
fillin $varname4tx $varname4confounder
predict E_outcome, pr
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
     +-------------------------------------------------------------------------------------+
     | qsmk   age_gt50     n   _fillin   E_outc~e   total   n_by_c~r   E_outc~x       diff |
     |-------------------------------------------------------------------------------------|
  1. |    0          0   868         0   .0823889    1629       1121   .1920047   .0108499 |
  2. |    0          1   333         0   .4338933    1629        508   .1920047   .0108499 |
  3. |    1          0   253         0   .0888793    1629       1121   .2028546   .0108499 |
  4. |    1          1   175         0   .4543631    1629        508   .2028546   .0108499 |
     +-------------------------------------------------------------------------------------+
*/

global varname4id seqn
global varname4outcome death
global varname4tx qsmk
global varname4confounder age_gt50 sex race exercise
global i_varname4confounder
foreach varname of global varname4confounder {
    global i_varname4confounder $i_varname4confounder i.`varname'
}
logit $varname4outcome $varname4tx $i_varname4confounder
preserve
collapse (count) n = $varname4id, by($varname4tx $varname4confounder)
fillin $varname4tx $varname4confounder
predict E_outcome, pr
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
     | qsmk   sex   race   exercise   age_gt50     n   _fillin   E_outc~e   total   n_by_c~r   E_outc~x       diff |
     |-------------------------------------------------------------------------------------------------------------|
  1. |    0     0      0          0          0   104         0   .0902793    1629        134   .1937839   .0047933 |
  2. |    0     0      0          1          0   143         0   .0917878    1629        204   .1937839   .0047933 |
  3. |    0     0      0          2          0    95         0   .1324816    1629        131   .1937839   .0047933 |
  4. |    0     0      1          0          0    11         0   .1172246    1629         12   .1937839   .0047933 |
  5. |    0     0      1          1          0    16         0   .1191244    1629         19   .1937839   .0047933 |
     |-------------------------------------------------------------------------------------------------------------|
  6. |    0     0      1          2          0    24         0   .1696736    1629         29   .1937839   .0047933 |
  7. |    0     1      0          0          0    79         0   .0484587    1629         92   .1937839   .0047933 |
  8. |    0     1      0          1          0   175         0   .0493063    1629        223   .1937839   .0047933 |
  9. |    0     1      0          2          0   144         0   .0726731    1629        185   .1937839   .0047933 |
 10. |    0     1      1          0          0     3         0   .0637973    1629          4   .1937839   .0047933 |
     |-------------------------------------------------------------------------------------------------------------|
 11. |    0     1      1          1          0    26         0   .0648949    1629         34   .1937839   .0047933 |
 12. |    0     1      1          2          0    48         0   .0949118    1629         54   .1937839   .0047933 |
 13. |    0     0      0          0          1    37         0   .4525754    1629         50   .1937839   .0047933 |
 14. |    0     0      0          1          1    56         0   .4570959    1629         92   .1937839   .0047933 |
 15. |    0     0      0          2          1    50         0   .5599024    1629         94   .1937839   .0047933 |
     |-------------------------------------------------------------------------------------------------------------|
 16. |    0     0      1          0          1     5         0   .5252236    1629          6   .1937839   .0047933 |
 17. |    0     0      1          1          1    10         0   .5297676    1629         13   .1937839   .0047933 |
 18. |    0     0      1          2          1    11         0   .6299534    1629         15   .1937839   .0047933 |
 19. |    0     1      0          0          1     8         0   .2978802    1629         19   .1937839   .0047933 |
 20. |    0     1      0          1          1    66         0   .3017071    1629         88   .1937839   .0047933 |
     |-------------------------------------------------------------------------------------------------------------|
 21. |    0     1      0          2          1    67         0   .3949921    1629        102   .1937839   .0047933 |
 22. |    0     1      1          0          1     .         1   .3621229    1629          0   .1937839   .0047933 |
 23. |    0     1      1          1          1     4         0   .3663446    1629          4   .1937839   .0047933 |
 24. |    0     1      1          2          1    19         0     .46627    1629         25   .1937839   .0047933 |
 25. |    1     0      0          0          0    30         0   .0934301    1629        134   .1985772   .0047933 |
     |-------------------------------------------------------------------------------------------------------------|
 26. |    1     0      0          1          0    61         0   .0949858    1629        204   .1985772   .0047933 |
 27. |    1     0      0          2          0    36         0   .1368837    1629        131   .1985772   .0047933 |
 28. |    1     0      1          0          0     1         0   .1211906    1629         12   .1985772   .0047933 |
 29. |    1     0      1          1          0     3         0   .1231457    1629         19   .1985772   .0047933 |
 30. |    1     0      1          2          0     5         0   .1750622    1629         29   .1985772   .0047933 |
     |-------------------------------------------------------------------------------------------------------------|
 31. |    1     1      0          0          0    13         0   .0502305    1629         92   .1985772   .0047933 |
 32. |    1     1      0          1          0    48         0   .0511074    1629        223   .1985772   .0047933 |
 33. |    1     1      0          2          0    41         0   .0752603    1629        185   .1985772   .0047933 |
 34. |    1     1      1          0          0     1         0   .0660911    1629          4   .1985772   .0047933 |
 35. |    1     1      1          1          0     8         0   .0672253    1629         34   .1985772   .0047933 |
     |-------------------------------------------------------------------------------------------------------------|
 36. |    1     1      1          2          0     6         0   .0982069    1629         54   .1985772   .0047933 |
 37. |    1     0      0          0          1    13         0   .4619499    1629         50   .1985772   .0047933 |
 38. |    1     0      0          1          1    36         0   .4664842    1629         92   .1985772   .0047933 |
 39. |    1     0      0          2          1    44         0   .5691886    1629         94   .1985772   .0047933 |
 40. |    1     0      1          0          1     1         0   .5346333    1629          6   .1985772   .0047933 |
     |-------------------------------------------------------------------------------------------------------------|
 41. |    1     0      1          1          1     3         0   .5391662    1629         13   .1985772   .0047933 |
 42. |    1     0      1          2          1     4         0   .6387153    1629         15   .1985772   .0047933 |
 43. |    1     1      0          0          1    11         0   .3058407    1629         19   .1985772   .0047933 |
 44. |    1     1      0          1          1    22         0   .3097247    1629         88   .1985772   .0047933 |
 45. |    1     1      0          2          1    35         0   .4040542    1629        102   .1985772   .0047933 |
     |-------------------------------------------------------------------------------------------------------------|
 46. |    1     1      1          0          1     .         1   .3708933    1629          0   .1985772   .0047933 |
 47. |    1     1      1          1          1     .         1   .3751571    1629          4   .1985772   .0047933 |
 48. |    1     1      1          2          1     6         0   .4756818    1629         25   .1985772   .0047933 |
     +-------------------------------------------------------------------------------------------------------------+
*/

