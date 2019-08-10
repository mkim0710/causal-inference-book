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


