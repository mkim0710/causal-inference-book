# function.data.ipw.PersonTime.SWsvyglmOutcome_Exposure_k.source.r
# https://github.com/mkim0710/causal-inference-book/blob/master/data/function.data.ipw.PersonTime.SWsvyglmOutcome_Exposure_k.source.r

library(tidyverse)
# nhefs = readRDS(gzcon(url("https://github.com/mkim0710/causal-inference-book/raw/master/data/nhefs.rds")))
# nhefs = readRDS(url("https://github.com/mkim0710/causal-inference-book/raw/master/data/nhefs.rds"))
# nhefs.surv = readRDS(url("https://github.com/mkim0710/causal-inference-book/raw/master/data/nhefs.surv.rds"))
nhefs.ipw.PersonTime = readRDS(url("https://github.com/mkim0710/causal-inference-book/raw/master/data/nhefs.ipw.PersonTime.rds"))



nhefs %>% select(seqn, death, matches("event"), matches("time")) %>% filter(seqn %in% c(21102, 14056, 2721)) #----
# nhefs.surv %>% select(seqn, death, matches("event"), matches("time")) %>% filter(seqn %in% c(21102, 14056, 2721)) #----
nhefs.ipw.PersonTime %>% select(seqn, survtime, death, PeriodSeq, Period, event, k, Dk_plus1) %>% filter(seqn %in% c(21102, 14056, 2721)) #----
# > nhefs %>% select(seqn, death, matches("event"), matches("time")) %>% filter(seqn %in% c(21102, 14056, 2721)) #----
# # A tibble: 3 x 3
#    seqn death survtime
#   <dbl> <dbl>    <dbl>
# 1  2721     1        2
# 2 14056     1        4
# 3 21102     1        1
# > nhefs.ipw.PersonTime %>% select(seqn, survtime, death, PeriodSeq, Period, event, k, Dk_plus1) %>% filter(seqn %in% c(21102, 14056, 2721)) #----
# # A tibble: 7 x 8
#    seqn survtime death PeriodSeq Period event     k Dk_plus1
#   <dbl>    <dbl> <dbl>     <int> <fct>  <lgl> <dbl> <lgl>   
# 1  2721        2     1         1 (0,1]  FALSE     0 FALSE   
# 2  2721        2     1         2 (1,2]  TRUE      1 TRUE    
# 3 14056        4     1         1 (0,1]  FALSE     0 FALSE   
# 4 14056        4     1         2 (1,2]  FALSE     1 FALSE   
# 5 14056        4     1         3 (2,3]  FALSE     2 FALSE   
# 6 14056        4     1         4 (3,4]  TRUE      3 TRUE    
# 7 21102        1     1         1 (0,1]  TRUE      0 TRUE    










#@@@@ Outcome Model with IP Weights -----
#@@@ MH) fit of weighted hazards model -----

#@ nhefs.ipw.PersonTime.SWsvyglmOutcome_Exposure_k ========
data = nhefs.ipw.PersonTime %>% select(Dk_plus1, Exposure, k, StabilizedWeight)
library(survey)
data.svydesign = data %>% svydesign(~1, weights = ~StabilizedWeight, data = .)
nhefs.ipw.PersonTime.SWsvyglmOutcome_Exposure_k =
    svyglm(formula = Dk_plus1 ~ Exposure * (k + I(k^2)), design = data.svydesign)
nhefs.ipw.PersonTime.SWsvyglmOutcome_Exposure_k %>% {cbind( `exp(coef(.))` = exp(coef(.)), exp(confint.default(.)), `Pr(>|t|)` = summary(.)$coefficients[,"Pr(>|t|)"] )} %>% round(2) %>% as.data.frame %>% rownames_to_column %>% as.tibble #----
# > nhefs.ipw.PersonTime.SWsvyglmOutcome_Exposure_k %>% {cbind( `exp(coef(.))` = exp(coef(.)), exp(confint.default(.)), `Pr(>|t|)` = summary(.)$coefficients[,"Pr(>|t|)"] )} %>% round(2) %>% as.data.frame %>% rownames_to_column %>% as.tibble #----
# # A tibble: 6 x 5
#   rowname         `exp(coef(.))` `2.5 %` `97.5 %` `Pr(>|t|)`
#   <chr>                    <dbl>   <dbl>    <dbl>      <dbl>
# 1 (Intercept)                  1       1        1       0   
# 2 Exposure                     1       1        1       0.64
# 3 k                            1       1        1       0.01
# 4 I(k^2)                       1       1        1       0.08
# 5 Exposure:k                   1       1        1       0.19
# 6 Exposure:I(k^2)              1       1        1       0.08





#@ MH) creation of dataset with all time points under each treatment level =====
data.PersonTime.glmOutcome_Exposure_k = nhefs.ipw.PersonTime.SWsvyglmOutcome_Exposure_k
data = nhefs.ipw.PersonTime %>% mutate(Exposure = qsmk, Dk_plus1 = event, k = k) 
data %>% select(k) %>% distinct %>% arrange(k) %>% {rbind(mutate(., Exposure = 0), mutate(., Exposure = 1))} %>% 
    mutate(pNoEvent_k = 1 - predict(data.PersonTime.glmOutcome_Exposure_k, newdata = ., type = "response")) %>% 
    group_by(Exposure) %>% mutate(pNoEvent_k.cumprod = pNoEvent_k %>% cumprod) %>% 
    rownames_to_column %>% group_by(Exposure) %>% do(tail(.,6)) %>% #----
as.tibble
# > nhefs.ipw.PersonTime %>% mutate(Exposure = qsmk, Dk_plus1 = event, k = k) %>% 
# +     select(k) %>% distinct %>% arrange(k) %>% 
# +     {rbind(mutate(., Exposure = 0), mutate(., Exposure = 1))} %>% 
# +     mutate(pNoEvent_k = 1 - predict(data.PersonTime.glmOutcome_Exposure_k, newdata = ., type = "response")) %>% 
# +     group_by(Exposure) %>% mutate(pNoEvent_k.cumprod = pNoEvent_k %>% cumprod) %>% 
# +     rownames_to_column %>% group_by(Exposure) %>% do(tail(.,6)) %>% #----
# +     as.tibble
# # A tibble: 12 x 5
#    rowname     k Exposure pNoEvent_k pNoEvent_k.cumprod
#    <chr>   <dbl>    <dbl>      <dbl>              <dbl>
#  1 115       114        0      0.998              0.812
#  2 116       115        0      0.998              0.811
#  3 117       116        0      0.998              0.809
#  4 118       117        0      0.998              0.808
#  5 119       118        0      0.998              0.806
#  6 120       119        0      0.998              0.805
#  7 235       114        1      0.999              0.810
#  8 236       115        1      0.999              0.809
#  9 237       116        1      0.999              0.809
# 10 238       117        1      0.999              0.808
# 11 239       118        1      0.999              0.807
# 12 240       119        1      0.999              0.807



#@ ggplot -----
data.PersonTime.glmOutcome_Exposure_k = nhefs.ipw.PersonTime.SWsvyglmOutcome_Exposure_k
data = nhefs.ipw.PersonTime %>% mutate(Exposure = qsmk, Dk_plus1 = event, k = k) 
g = data %>% select(k) %>% distinct %>% arrange(k) %>% {rbind(mutate(., Exposure = 0), mutate(., Exposure = 1))} %>% 
    mutate(pNoEvent_k = 1 - predict(data.PersonTime.glmOutcome_Exposure_k, newdata = ., type = "response")) %>% 
    group_by(Exposure) %>% mutate(pNoEvent_k.cumprod = pNoEvent_k %>% cumprod) %>% 
    ungroup %>% mutate(Exposure = Exposure %>% as.factor) %>% 
    ggplot(aes(x = k, y = pNoEvent_k.cumprod, linetype = Exposure, group = Exposure)) + 
    geom_line()


g
g+theme_classic()
g+theme_bw()
g+theme_minimal()+theme(legend.position="bottom")
g+theme_bw()+theme(legend.position="bottom")

filename = paste0("nhefs.ipw.PersonTime.SWsvyglmOutcome_Exposure_k", ".ggplot")
g+theme_bw()+theme(legend.position="bottom")
ggsave(paste0(filename, ".pdf"), width = 8, height = 6)
ggsave(paste0(filename, ".png"), width = 8, height = 6)



#@ (cumulative incidence) -----
data.PersonTime.glmOutcome_Exposure_k = nhefs.ipw.PersonTime.SWsvyglmOutcome_Exposure_k
data = nhefs.ipw.PersonTime %>% mutate(Exposure = qsmk, Dk_plus1 = event, k = k) 
g = data %>% select(k) %>% distinct %>% arrange(k) %>% {rbind(mutate(., Exposure = 0), mutate(., Exposure = 1))} %>% 
    mutate(pNoEvent_k = 1 - predict(data.PersonTime.glmOutcome_Exposure_k, newdata = ., type = "response")) %>% 
    group_by(Exposure) %>% mutate(pNoEvent_k.cumprod = pNoEvent_k %>% cumprod) %>% 
    ungroup %>% mutate(Exposure = Exposure %>% as.factor) %>% 
    ggplot(aes(x = k, y = 1 - pNoEvent_k.cumprod, linetype = Exposure, group = Exposure)) + 
    geom_line()

filename = paste0("nhefs.ipw.PersonTime.SWsvyglmOutcome_Exposure_k", ".ggplot", " (cumulative incidence)")
g+theme_bw()+theme(legend.position="bottom")
ggsave(paste0(filename, ".pdf"), width = 8, height = 6)
ggsave(paste0(filename, ".png"), width = 8, height = 6)




#@ end ----
