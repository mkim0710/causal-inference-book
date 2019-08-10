nhefs = readRDS(gzcon(url("https://github.com/mkim0710/causal-inference-book/raw/master/data/nhefs.rds")))
nhefs = readRDS(url("https://github.com/mkim0710/causal-inference-book/raw/master/data/nhefs.rds"))

nhefs.surv = readRDS(url("https://github.com/mkim0710/causal-inference-book/raw/master/data/nhefs.surv.rds"))



nhefs %>% select(seqn, death, matches("event"), matches("time")) %>% filter(seqn %in% c(21102, 14056, 2721)) #----
nhefs.surv %>% select(seqn, death, matches("event"), matches("time")) %>% filter(seqn %in% c(21102, 14056, 2721)) #----
# > nhefs %>% select(seqn, death, matches("event"), matches("time")) %>% filter(seqn %in% c(21102, 14056, 2721)) #----
# # A tibble: 3 x 3
#    seqn death survtime
#   <dbl> <dbl>    <dbl>
# 1  2721     1        2
# 2 14056     1        4
# 3 21102     1        1
# > nhefs.surv %>% select(seqn, death, matches("event"), matches("time")) %>% filter(seqn %in% c(21102, 14056, 2721)) #----
# # A tibble: 7 x 6
#    seqn death event survtime  time timesq
#   <dbl> <dbl> <dbl>    <dbl> <dbl>  <dbl>
# 1  2721     1     0        2     0      0
# 2  2721     1     1        2     1      1
# 3 14056     1     0        4     0      0
# 4 14056     1     0        4     1      1
# 5 14056     1     0        4     2      4
# 6 14056     1     1        4     3      9
# 7 21102     1     1        1     0      0


all.equal(
    nhefs.surv
    ,
    nhefs %>% mutate(time = survtime %>% map(function(x) 1:x - 1)) %>% unnest %>%
        mutate(
            event = as.numeric(time==survtime-1 & death==1)
            , timesq = time*time
        )
)
# > all.equal(
# +     nhefs.surv
# +     ,
# +     nhefs %>% mutate(time = survtime %>% map(function(x) 1:x - 1)) %>% unnest %>%
# +         mutate(
# +             event = as.numeric(time==survtime-1 & death==1)
# +             , timesq = time*time
# +         )
# + )
# [1] TRUE


Interval = 50
nhefs %>%
    mutate(PeriodSeq = survtime %>% map(function(x) 1L:ceiling(x/Interval))) %>% unnest %>%  
    mutate(
        Period = paste0("(", (PeriodSeq-1)*Interval, ",", PeriodSeq*Interval, "]") %>% as.factor
        , time = PeriodSeq * Interval
        , timesq = time*time
        , event = (death == 1) & (survtime <= PeriodSeq * Interval)
        , k = PeriodSeq - 1  # defined as in hernanrobins_v2.17.22 $17.2 From hazards to risks
        , ksq = k*k
        , Dk_plus1 = event  # defined as in hernanrobins_v2.17.22 $17.2 From hazards to risks 
        , Exposure = qsmk
        , `Exposure:k` = Exposure * k
        , `I(Exposure*ksq)` = Exposure * ksq
    ) %>% 
    select(seqn, survtime, death, PeriodSeq, Period, time, timesq, event, k, Dk_plus1) #----
# > Interval = 50
# > nhefs %>%
# +     mutate(PeriodSeq = survtime %>% map(function(x) 1L:ceiling(x/Interval))) %>% unnest %>%  
# +     mutate(
# +         Period = paste0("(", (PeriodSeq-1)*Interval, ",", PeriodSeq*Interval, "]") %>% as.factor
# +         , time = PeriodSeq * Interval
# +         , event = (death == 1) & (survtime <= PeriodSeq * Interval)
# +         , timesq = time*time
# +         , k = PeriodSeq - 1  # defined as in hernanrobins_v2.17.22 $17.2 From hazards to risks
# +         , Dk_plus1 = event  # defined as in hernanrobins_v2.17.22 $17.2 From hazards to risks 
# +     ) %>% 
# +     select(seqn, survtime, death, PeriodSeq, Period, time, timesq, event, k, Dk_plus1) #----
# # A tibble: 4,480 x 10
#     seqn survtime death PeriodSeq Period     time timesq event     k Dk_plus1
#    <dbl>    <dbl> <dbl>     <int> <fct>     <dbl>  <dbl> <lgl> <dbl> <lgl>   
#  1   233      120     0         1 (0,50]       50   2500 FALSE     0 FALSE   
#  2   233      120     0         2 (50,100]    100  10000 FALSE     1 FALSE   
#  3   233      120     0         3 (100,150]   150  22500 FALSE     2 FALSE   
#  4   235      120     0         1 (0,50]       50   2500 FALSE     0 FALSE   
#  5   235      120     0         2 (50,100]    100  10000 FALSE     1 FALSE   
#  6   235      120     0         3 (100,150]   150  22500 FALSE     2 FALSE   
#  7   244      120     0         1 (0,50]       50   2500 FALSE     0 FALSE   
#  8   244      120     0         2 (50,100]    100  10000 FALSE     1 FALSE   
#  9   244      120     0         3 (100,150]   150  22500 FALSE     2 FALSE   
# 10   245       26     1         1 (0,50]       50   2500 TRUE      0 TRUE    
# # ... with 4,470 more rows











#@@@ MH) coxph() =====

#@ nhefs.surv.glm_Exposure ====
data = nhefs.surv %>% mutate(Exposure = qsmk, Dk_plus1 = event, k = time) %>% 
  select(
    Dk_plus1, Exposure, k
  )
nhefs.surv.glm_Exposure = glm(formula = Dk_plus1 ~ Exposure + (k + I(k^2)) + . , data = data, family = binomial)
nhefs.surv.glm_Exposure %>% {cbind( coef(.), confint.default(.) )} %>% exp %>% round(2) #----
# nhefs.surv.glm_Exposure %>% summary #----
# > nhefs.surv.glm_Exposure %>% {cbind( coef(.), confint.default(.) )} %>% exp %>% round(2) #----
#                  2.5 % 97.5 %
# (Intercept) 0.00  0.00   0.00
# Exposure    1.39  1.10   1.77
# k           1.02  1.01   1.04
# I(k^2)      1.00  1.00   1.00


#@ nhefs.coxph_Exposure ====
data = nhefs %>% mutate(Exposure = qsmk, event = death, time = survtime) %>% 
  select(
    time, event, Exposure
  )
library(survival)
nhefs.coxph_Exposure = coxph(formula = Surv(time = time, event = event) ~ . , data = data, method = "breslow")
nhefs.coxph_Exposure %>% {cbind( coef(.), confint(.) )} %>% exp %>% round(2) #----
# nhefs.coxph_Exposure %>% summary #----
# > nhefs.coxph_Exposure %>% {cbind( coef(.), confint(.) )} %>% exp %>% round(2) #----
#               2.5 % 97.5 %
# Exposure 1.39   1.1   1.76






                                      




#@ -----
#@ MH) fit of parametric hazards model ========
data = nhefs.surv %>% mutate(Exposure = qsmk, Dk_plus1 = event, k = time) %>% 
  select(
    Dk_plus1, Exposure, k
  )
nhefs.surv.glm_Exposure_k = glm(formula = Dk_plus1 ~ Exposure * (k + I(k^2)) + . , data = data, family = binomial)
nhefs.surv.glm_Exposure_k %>% {cbind( coef(.), confint.default(.) )} %>% exp %>% round(2) #----
nhefs.surv.glm_Exposure_k %>% summary #----
# > nhefs.surv.glm_Exposure_k %>% {cbind( coef(.), confint.default(.) )} %>% exp %>% round(2) #----
#                      2.5 % 97.5 %
# (Intercept)     0.00  0.00   0.00
# Exposure        1.40  0.64   3.05
# k               1.02  1.00   1.04
# I(k^2)          1.00  1.00   1.00
# Exposure:k      1.01  0.98   1.04
# Exposure:I(k^2) 1.00  1.00   1.00


#@ MH) creation of dataset with all time points under each treatment level =====
data.PersonTime.glm_Exposure_k = nhefs.surv.glm_Exposure_k
nhefs.surv %>% mutate(Exposure = qsmk, Dk_plus1 = event, k = time) %>% 
  select(k) %>% distinct %>% arrange(k) %>% 
  {rbind(mutate(., Exposure = 0), mutate(., Exposure = 1))} %>% 
  mutate(pNoEvent_k = 1 - predict(data.PersonTime.glm_Exposure_k, newdata = ., type = "response")) %>% 
  group_by(Exposure) %>% mutate(pNoEvent_k.cumprod = pNoEvent_k %>% cumprod) %>% 
  as.tibble
# > nhefs.surv %>% mutate(Exposure = qsmk, Dk_plus1 = event, k = time) %>% 
# +   select(k) %>% distinct %>% arrange(k) %>% 
# +   {rbind(mutate(., Exposure = 0), mutate(., Exposure = 1))} %>% 
# +   mutate(pNoEvent_k = 1 - predict(nhefs.surv.glm_Exposure_k, newdata = ., type = "response")) %>% 
# +   group_by(Exposure) %>% mutate(pNoEvent_k.cumprod = pNoEvent_k %>% cumprod) %>% 
# +   as.tibble
# # A tibble: 240 x 4
#        k Exposure pNoEvent_k pNoEvent_k.cumprod
#    <dbl>    <dbl>      <dbl>              <dbl>
#  1     0        0      0.999              0.999
#  2     1        0      0.999              0.998
#  3     2        0      0.999              0.997
#  4     3        0      0.999              0.996
#  5     4        0      0.999              0.995
#  6     5        0      0.999              0.994
#  7     6        0      0.999              0.993
#  8     7        0      0.999              0.992
#  9     8        0      0.999              0.991
# 10     9        0      0.999              0.990
# # ... with 230 more rows


data.PersonTime.glm_Exposure_k = nhefs.surv.glm_Exposure_k
nhefs.surv %>% mutate(Exposure = qsmk, Dk_plus1 = event, k = time) %>% 
  select(k) %>% distinct %>% arrange(k) %>% 
  {rbind(mutate(., Exposure = 0), mutate(., Exposure = 1))} %>% 
  mutate(pNoEvent_k = 1 - predict(data.PersonTime.glm_Exposure_k, newdata = ., type = "response")) %>% 
  group_by(Exposure) %>% mutate(pNoEvent_k.cumprod = pNoEvent_k %>% cumprod) %>% 
  ungroup %>% mutate(Exposure = Exposure %>% as.factor) %>% 
  ggplot(aes(x = k, y = pNoEvent_k.cumprod, color = Exposure, group = Exposure)) + 
  geom_line()



