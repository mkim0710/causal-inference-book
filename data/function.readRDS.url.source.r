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
        , event = (death == 1) & (survtime <= PeriodSeq * Interval)
        , timesq = time*time
        , k = PeriodSeq - 1  # defined as in hernanrobins_v2.17.22 $17.2 From hazards to risks
        , Dk_plus1 = event  # defined as in hernanrobins_v2.17.22 $17.2 From hazards to risks 
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
