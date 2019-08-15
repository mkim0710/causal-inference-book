# nhefs.ipw.PersonTime from nhefs.r
# https://github.com/mkim0710/causal-inference-book/blob/master/data/nhefs.ipw.PersonTime from nhefs.r


# nhefs = readRDS(gzcon(url("https://github.com/mkim0710/causal-inference-book/raw/master/data/nhefs.rds")))
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














#@@@@ Exposure Propensity Model -----

# Unstabilized Weight: 1/f(A|L) = 1 / { A·P(A=1|L) + (1-A)·(1-P(A=1|L)) }
# Stabilized Weight: f(A)/f(A|L) = { A·P(A=1) + (1-A)·(1-P(A=1)) } / { A·P(A=1|L) + (1-A)·(1-P(A=1|L)) }
# Stabilized Weight: f(A)/f(A|L) = A·{P(A=1)/P(A=1|L)} + (1-A)·{(1-P(A=1))/(1-P(A=1|L))}

#@@ MH) estimation of denominator of ip weights - this is enough for unstabilized weights: : 1/f(A|L) = 1 / { A·P(A=1|L) + (1-A)·(1-P(A=1|L)) } ----
nhefs.glmExposure_Covariate = 
    glm(qsmk ~ sex + race + age + I(age*age) + as.factor(education)
        + smokeintensity + I(smokeintensity*smokeintensity)
        + smokeyrs + I(smokeyrs*smokeyrs) + as.factor(exercise)
        + as.factor(active) + wt71 + I(wt71*wt71), 
        data=nhefs, family=binomial())
nhefs.glmExposure_Covariate %>% {cbind( `exp(coef(.))` = exp(coef(.)), exp(confint.default(.)), `Pr(>|z|)` = summary(.)$coefficients[,"Pr(>|z|)"] )} %>% round(2) %>% as.data.frame %>% rownames_to_column %>% as.tibble #----
# > nhefs.glmExposure_Covariate %>% {cbind( `exp(coef(.))` = exp(coef(.)), exp(confint.default(.)), `Pr(>|z|)` = summary(.)$coefficients[,"Pr(>|z|)"] )} %>% round(2) %>% as.data.frame %>% rownames_to_column %>% as.tibble #----
# # A tibble: 19 x 5
#    rowname                            `exp(coef(.))` `2.5 %` `97.5 %` `Pr(>|z|)`
#    <chr>                                       <dbl>   <dbl>    <dbl>      <dbl>
#  1 (Intercept)                                  0.14   0.01      1.56       0.11
#  2 sex                                          0.6    0.45      0.8        0   
#  3 race                                         0.43   0.290     0.64       0   
#  4 age                                          1.11   1.01      1.22       0.04
#  5 I(age * age)                                 1      1         1          0.23
#  6 as.factor(education)2                        0.91   0.62      1.32       0.61
#  7 as.factor(education)3                        1.02   0.73      1.42       0.93
#  8 as.factor(education)4                        0.96   0.570     1.61       0.87
#  9 as.factor(education)5                        1.46   0.95      2.25       0.08
# 10 smokeintensity                               0.94   0.91      0.96       0   
# 11 I(smokeintensity * smokeintensity)           1      1         1          0   
# 12 smokeyrs                                     0.93   0.88      0.98       0.01
# 13 I(smokeyrs * smokeyrs)                       1      1         1          0.06
# 14 as.factor(exercise)1                         1.34   0.95      1.88       0.09
# 15 as.factor(exercise)2                         1.43   1         2.03       0.05
# 16 as.factor(active)1                           1.01   0.78      1.3        0.93
# 17 as.factor(active)2                           1.07   0.71      1.61       0.74
# 18 wt71                                         0.99   0.95      1.03       0.56
# 19 I(wt71 * wt71)                               1      1         1          0.37




#@@ MH) computation of estimated weights ----
# Unstabilized Weight: 1/f(A|L) = 1 / { A·P(A=1|L) + (1-A)·(1-P(A=1|L)) }
# Stabilized Weight: f(A)/f(A|L) = { A·P(A=1) + (1-A)·(1-P(A=1)) } / { A·P(A=1|L) + (1-A)·(1-P(A=1|L)) }
# Stabilized Weight: f(A)/f(A|L) = A·{P(A=1)/P(A=1|L)} + (1-A)·{(1-P(A=1))/(1-P(A=1|L))}

data.glmExposure_Covariate = nhefs.glmExposure_Covariate
nhefs %>% transmute(
    Exposure = qsmk
    , pExposure = mean(Exposure)
    , pExposure_Covariate = predict(data.glmExposure_Covariate, newdata = ., type = "response")
    , UnstabilizedWeight = if_else(Exposure==1, 1/pExposure_Covariate, 1/(1-pExposure_Covariate))
    , StabilizedWeight = if_else(Exposure==1, pExposure/pExposure_Covariate, (1-pExposure)/(1-pExposure_Covariate))
) %>% rownames_to_column %>% group_by(Exposure) %>% do(tail(.,3)) #----
# > nhefs %>% transmute(
# +     Exposure = qsmk
# +     , pExposure = mean(Exposure)
# +     , pExposure_Covariate = predict(data.glmExposure_Covariate, newdata = ., type = "response")
# +     , UnstabilizedWeight = if_else(Exposure==1, 1/pExposure_Covariate, 1/(1-pExposure_Covariate))
# +     , StabilizedWeight = if_else(Exposure==1, pExposure/pExposure_Covariate, (1-pExposure)/(1-pExposure_Covariate))
# + ) %>% rownames_to_column %>% group_by(Exposure) %>% do(tail(.,3)) #----
# # A tibble: 6 x 6
# # Groups:   Exposure [2]
#   rowname Exposure pExposure pExposure_Covariate UnstabilizedWeight StabilizedWeight
#   <chr>      <dbl>     <dbl>               <dbl>              <dbl>            <dbl>
# 1 1626           0     0.263               0.146               1.17            0.864
# 2 1627           0     0.263               0.138               1.16            0.855
# 3 1628           0     0.263               0.448               1.81            1.34 
# 4 1606           1     0.263               0.793               1.26            0.331
# 5 1607           1     0.263               0.239               4.19            1.10 
# 6 1629           1     0.263               0.178               5.63            1.48 


# data.glmExposure_Covariate = nhefs.glmExposure_Covariate
# nhefs =
#     nhefs %>% mutate(
#         Exposure = qsmk
#         , pExposure = mean(Exposure)
#         , pExposure_Covariate = predict(data.glmExposure_Covariate, newdata = ., type = "response")
#         , UnstabilizedWeight = if_else(Exposure==1, 1/pExposure_Covariate, 1/(1-pExposure_Covariate))
#         , StabilizedWeight = if_else(Exposure==1, pExposure/pExposure_Covariate, (1-pExposure)/(1-pExposure_Covariate))
#     )










#@@@ creation of person-month data -----
nhefs.ipw <- splitstackshape::expandRows(nhefs, "survtime", drop=F) 
nhefs.ipw$time <- sequence(rle(nhefs.ipw$seqn)$lengths)-1
nhefs.ipw$event <- ifelse(nhefs.ipw$time==nhefs.ipw$survtime-1 & 
                              nhefs.ipw$death==1, 1, 0)
nhefs.ipw$timesq <- nhefs.ipw$time^2

identical(
    nhefs.surv %>% select(seqn, death, yrdth, modth, survtime, time)
    ,
    nhefs.ipw %>% select(seqn, death, yrdth, modth, survtime, time)
)
# > identical(
# +     nhefs.surv %>% select(seqn, death, yrdth, modth, survtime, time)
# +     ,
# +     nhefs.ipw %>% select(seqn, death, yrdth, modth, survtime, time)
# + )
# [1] TRUE

#@ end ----
write_rds(nhefs.ipw, "nhefs.ipw.rds", "gz", compression = 9)





#@@@ MH) creation of person-month data -----
Interval = 1
nhefs %>%
    mutate(PeriodSeq = survtime %>% map(function(x) 1L:ceiling(x/Interval))) %>% unnest %>%  
    mutate(
        Period = paste0("(", (PeriodSeq-1)*Interval, ",", PeriodSeq*Interval, "]") %>% as.factor
        # , time = PeriodSeq * Interval
        # , timesq = time*time
        , event = (death == 1) & (survtime <= PeriodSeq * Interval)
        , k = PeriodSeq - 1  # defined as in hernanrobins_v2.17.22 $17.2 From hazards to risks
        , ksq = k*k
        , Dk_plus1 = event  # defined as in hernanrobins_v2.17.22 $17.2 From hazards to risks 
        , Exposure = qsmk
        , `Exposure:k` = Exposure * k
        , `I(Exposure*ksq)` = Exposure * ksq
    ) %>% 
    # select(seqn, survtime, death, PeriodSeq, Period, time, timesq, event, k, Dk_plus1) #----
select(seqn, survtime, death, PeriodSeq, Period, event, k, Dk_plus1) #----
# > nhefs %>%
# +     mutate(PeriodSeq = survtime %>% map(function(x) 1L:ceiling(x/Interval))) %>% unnest %>%  
# +     mutate(
# +         Period = paste0("(", (PeriodSeq-1)*Interval, ",", PeriodSeq*Interval, "]") %>% as.factor
# +         # , time = PeriodSeq * Interval
# +         # , timesq = time*time
# +         , event = (death == 1) & (survtime <= PeriodSeq * Interval)
# +         , k = PeriodSeq - 1  # defined as in hernanrobins_v2.17.22 $17.2 From hazards to risks
# +         , ksq = k*k
# +         , Dk_plus1 = event  # defined as in hernanrobins_v2.17.22 $17.2 From hazards to risks 
# +         , Exposure = qsmk
# +         , `Exposure:k` = Exposure * k
# +         , `I(Exposure*ksq)` = Exposure * ksq
# +     ) %>% 
# +     # select(seqn, survtime, death, PeriodSeq, Period, time, timesq, event, k, Dk_plus1) #----
# +     select(seqn, survtime, death, PeriodSeq, Period, event, k, Dk_plus1) #----
# # A tibble: 176,764 x 8
#     seqn survtime death PeriodSeq Period event     k Dk_plus1
#    <dbl>    <dbl> <dbl>     <int> <fct>  <lgl> <dbl> <lgl>   
#  1   233      120     0         1 (0,1]  FALSE     0 FALSE   
#  2   233      120     0         2 (1,2]  FALSE     1 FALSE   
#  3   233      120     0         3 (2,3]  FALSE     2 FALSE   
#  4   233      120     0         4 (3,4]  FALSE     3 FALSE   
#  5   233      120     0         5 (4,5]  FALSE     4 FALSE   
#  6   233      120     0         6 (5,6]  FALSE     5 FALSE   
#  7   233      120     0         7 (6,7]  FALSE     6 FALSE   
#  8   233      120     0         8 (7,8]  FALSE     7 FALSE   
#  9   233      120     0         9 (8,9]  FALSE     8 FALSE   
# 10   233      120     0        10 (9,10] FALSE     9 FALSE   
# # ... with 176,754 more rows

nhefs.surv %>% select(seqn, death, matches("event"), matches("time")) #----
# > nhefs.surv %>% select(seqn, death, matches("event"), matches("time")) #----
# # A tibble: 176,764 x 6
#     seqn death event survtime  time timesq
#    <dbl> <dbl> <dbl>    <dbl> <dbl>  <dbl>
#  1   233     0     0      120     0      0
#  2   233     0     0      120     1      1
#  3   233     0     0      120     2      4
#  4   233     0     0      120     3      9
#  5   233     0     0      120     4     16
#  6   233     0     0      120     5     25
#  7   233     0     0      120     6     36
#  8   233     0     0      120     7     49
#  9   233     0     0      120     8     64
# 10   233     0     0      120     9     81
# # ... with 176,754 more rows



all.equal(
    nhefs.surv %>% select(seqn, death, matches("event"), matches("time"))
    , 
    nhefs %>%
        mutate(PeriodSeq = survtime %>% map(function(x) 1L:ceiling(x/Interval))) %>% unnest %>%  
        mutate(
            Period = paste0("(", (PeriodSeq-1)*Interval, ",", PeriodSeq*Interval, "]") %>% as.factor
            # , time = PeriodSeq * Interval
            # , timesq = time*time
            , event = (death == 1) & (survtime <= PeriodSeq * Interval)
            , k = PeriodSeq - 1  # defined as in hernanrobins_v2.17.22 $17.2 From hazards to risks
            , ksq = k*k
            , Dk_plus1 = event  # defined as in hernanrobins_v2.17.22 $17.2 From hazards to risks 
            , Exposure = qsmk
            , `Exposure:k` = Exposure * k
            , `I(Exposure*ksq)` = Exposure * ksq
        ) %>% 
        # select(seqn, survtime, death, PeriodSeq, Period, time, timesq, event, k, Dk_plus1) #----
    # select(seqn, survtime, death, PeriodSeq, Period, event, k, Dk_plus1) #----
    transmute(
        seqn = seqn
        , death = death
        , event = event %>% as.numeric
        , survtime = survtime
        , time = k
        , timesq = k*k
    )
)
# [1] TRUE






#@@ MH) PersonTime data transformation with weights ----
# Unstabilized Weight: 1/f(A|L) = 1 / { A·P(A=1|L) + (1-A)·(1-P(A=1|L)) }
# Stabilized Weight: f(A)/f(A|L) = { A·P(A=1) + (1-A)·(1-P(A=1)) } / { A·P(A=1|L) + (1-A)·(1-P(A=1|L)) }
# Stabilized Weight: f(A)/f(A|L) = A·{P(A=1)/P(A=1|L)} + (1-A)·{(1-P(A=1))/(1-P(A=1|L))}

#@ nhefs.ipw.PersonTime =====
Interval = 1
data.glmExposure_Covariate = nhefs.glmExposure_Covariate
nhefs.ipw.PersonTime =
    nhefs %>% mutate(
        Exposure = qsmk
        , pExposure = mean(Exposure)
        , pExposure_Covariate = predict(data.glmExposure_Covariate, newdata = ., type = "response")
        , UnstabilizedWeight = if_else(Exposure==1, 1/pExposure_Covariate, 1/(1-pExposure_Covariate))
        , StabilizedWeight = if_else(Exposure==1, pExposure/pExposure_Covariate, (1-pExposure)/(1-pExposure_Covariate))
    ) %>%
    mutate(PeriodSeq = survtime %>% map(function(x) 1L:ceiling(x/Interval))) %>% unnest %>%  
    mutate(
        Period = paste0("(", (PeriodSeq-1)*Interval, ",", PeriodSeq*Interval, "]") %>% as.factor
        # , time = PeriodSeq * Interval
        # , timesq = time*time
        , event = (death == 1) & (survtime <= PeriodSeq * Interval)
        , k = PeriodSeq - 1  # defined as in hernanrobins_v2.17.22 $17.2 From hazards to risks
        , ksq = k*k
        , Dk_plus1 = event  # defined as in hernanrobins_v2.17.22 $17.2 From hazards to risks 
        , Exposure = qsmk
        , `Exposure:k` = Exposure * k
        , `I(Exposure*ksq)` = Exposure * ksq
    )

nhefs.ipw.PersonTime %>% str #----
# > nhefs.ipw.PersonTime %>% str #----
# Classes ‘tbl_df’, ‘tbl’ and 'data.frame':	176764 obs. of  81 variables:
#  $ seqn               : num  233 233 233 233 233 233 233 233 233 233 ...
#  $ qsmk               : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ death              : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ yrdth              : num  NA NA NA NA NA NA NA NA NA NA ...
#  $ modth              : num  NA NA NA NA NA NA NA NA NA NA ...
#  $ dadth              : num  NA NA NA NA NA NA NA NA NA NA ...
#  $ sbp                : num  175 175 175 175 175 175 175 175 175 175 ...
#  $ dbp                : num  96 96 96 96 96 96 96 96 96 96 ...
#  $ sex                : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ age                : num  42 42 42 42 42 42 42 42 42 42 ...
#  $ race               : num  1 1 1 1 1 1 1 1 1 1 ...
#  $ income             : num  19 19 19 19 19 19 19 19 19 19 ...
#  $ marital            : num  2 2 2 2 2 2 2 2 2 2 ...
#  $ school             : num  7 7 7 7 7 7 7 7 7 7 ...
#  $ education          : num  1 1 1 1 1 1 1 1 1 1 ...
#  $ ht                 : num  174 174 174 174 174 ...
#  $ wt71               : num  79 79 79 79 79 ...
#  $ wt82               : num  68.9 68.9 68.9 68.9 68.9 ...
#  $ wt82_71            : num  -10.1 -10.1 -10.1 -10.1 -10.1 ...
#  $ birthplace         : num  47 47 47 47 47 47 47 47 47 47 ...
#  $ smokeintensity     : num  30 30 30 30 30 30 30 30 30 30 ...
#  $ smkintensity82_71  : num  -10 -10 -10 -10 -10 -10 -10 -10 -10 -10 ...
#  $ smokeyrs           : num  29 29 29 29 29 29 29 29 29 29 ...
#  $ asthma             : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ bronch             : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ tb                 : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ hf                 : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ hbp                : num  1 1 1 1 1 1 1 1 1 1 ...
#  $ pepticulcer        : num  1 1 1 1 1 1 1 1 1 1 ...
#  $ colitis            : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ hepatitis          : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ chroniccough       : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ hayfever           : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ diabetes           : num  1 1 1 1 1 1 1 1 1 1 ...
#  $ polio              : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ tumor              : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ nervousbreak       : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ alcoholpy          : num  1 1 1 1 1 1 1 1 1 1 ...
#  $ alcoholfreq        : num  1 1 1 1 1 1 1 1 1 1 ...
#  $ alcoholtype        : num  3 3 3 3 3 3 3 3 3 3 ...
#  $ alcoholhowmuch     : num  7 7 7 7 7 7 7 7 7 7 ...
#  $ pica               : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ headache           : num  1 1 1 1 1 1 1 1 1 1 ...
#  $ otherpain          : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ weakheart          : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ allergies          : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ nerves             : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ lackpep            : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ hbpmed             : num  1 1 1 1 1 1 1 1 1 1 ...
#  $ boweltrouble       : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ wtloss             : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ infection          : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ active             : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ exercise           : num  2 2 2 2 2 2 2 2 2 2 ...
#  $ birthcontrol       : num  2 2 2 2 2 2 2 2 2 2 ...
#  $ pregnancies        : num  NA NA NA NA NA NA NA NA NA NA ...
#  $ cholesterol        : num  197 197 197 197 197 197 197 197 197 197 ...
#  $ hightax82          : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ price71            : num  2.18 2.18 2.18 2.18 2.18 ...
#  $ price82            : num  1.74 1.74 1.74 1.74 1.74 ...
#  $ tax71              : num  1.1 1.1 1.1 1.1 1.1 ...
#  $ tax82              : num  0.462 0.462 0.462 0.462 0.462 ...
#  $ price71_82         : num  0.444 0.444 0.444 0.444 0.444 ...
#  $ tax71_82           : num  0.64 0.64 0.64 0.64 0.64 ...
#  $ survtime           : num  120 120 120 120 120 120 120 120 120 120 ...
#  $ pd.qsmk            : num  0.109 0.109 0.109 0.109 0.109 ...
#  $ pExposure_Covariate: num  0.109 0.109 0.109 0.109 0.109 ...
#  $ pn.qsmk            : num  0.263 0.263 0.263 0.263 0.263 ...
#  $ sw.a               : num  0.827 0.827 0.827 0.827 0.827 ...
#  $ Exposure           : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ pExposure          : num  0.263 0.263 0.263 0.263 0.263 ...
#  $ UnstabilizedWeight : num  1.12 1.12 1.12 1.12 1.12 ...
#  $ StabilizedWeight   : num  0.827 0.827 0.827 0.827 0.827 ...
#  $ PeriodSeq          : int  1 2 3 4 5 6 7 8 9 10 ...
#  $ Period             : Factor w/ 120 levels "(0,1]","(1,2]",..: 1 2 33 44 55 66 77 88 99 110 ...
#  $ event              : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
#  $ k                  : num  0 1 2 3 4 5 6 7 8 9 ...
#  $ ksq                : num  0 1 4 9 16 25 36 49 64 81 ...
#  $ Dk_plus1           : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
#  $ Exposure:k         : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ I(Exposure*ksq)    : num  0 0 0 0 0 0 0 0 0 0 ...



#@ end ----
write_rds(nhefs.ipw.PersonTime, "nhefs.ipw.PersonTime.rds", "gz", compression = 9)
