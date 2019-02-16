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
regress $varname4outcome $varname4tx $varname4confounder
restore, preserve
preserve
collapse (count) n = $varname4id, by($varname4tx $varname4confounder)
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


