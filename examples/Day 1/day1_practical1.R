#################################################
### Bayesian causal inference day 1 practical ###
#################################################

#################
### Section 1 ###
#################

set.seed(1)
setwd("~/Downloads/LeuvenBCI-main") # change filepath accordingly
load("data/lalonde.RData")

## extract variables
y <- ldw[,"re78"]
t <- ldw[,"treat"]
