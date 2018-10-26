library(data.table)

tmp = fread('./patient-therapist.csv')
str(tmp)

library(depmixS4)

norm.HMM.Generate <- function(n, m, mu, sd, gamma1, delta=NULL){
  if(m != length(gamma1[,1]) | m != length(gamma1[1,])){
    stop("The number of states given by m does not match the number of rows
and/or colums of the transition probability matrix gamma")
  }
  state <- numeric(n)
  # first attribute the first observation in the sequence
  if (is.null(delta)) {
    state[1] <- sample(1:m, 1, prob = solve(t(diag(m) - gamma1 + 1), rep(1, m)))
  } else {
    state[1] <- sample(1:m, 1, prob = delta)
  }
  # conditional on the first observation, sample the following observations
  for (i in 2:n){
    state[i] <- sample(1:m, 1, prob = gamma1[state[i-1],])
  }
  x <- rnorm(n, mean = mu[state], sd = sd[state])
  return(list(state = state, seq.x = x))
}
sample.data <- norm.HMM.Generate(n = 1000, m = 2, mu = c(3, 10), sd = c(1.5, 3),
                                 gamma1 = matrix(c(.8, .2, .2, .8), byrow = T, ncol = 2))

plot(sample.data$seq.x)

library(depmixS4)

model <- depmix(response = seq.x ~ 1,
                data=as.data.frame(sample.data),
                nstates=2,
                family=gaussian(),
                ntimes=1000)

fitted_model <- fit(model)
fitted_model
# degrees of freedom 7
# 4 - 2 each for the distributions - mu1, mu2, sd1, sd2
# 2 - transition probabilities 4, but each has to sum to 1, so effectively only 2
# 1 - initial state
summary(fitted_model)

posterior(fitted_model)
# output 3 columns - 1) most likely state, probability of each state 2 & 3

# -----------------

earthq.data <- c(13, 14, 8, 10, 16, 26, 32, 27, 18, 32, 36, 24, 22, 23, 22, 18, 25,
                 21, 21, 14, 8, 11, 14, 23, 18, 17, 19, 20, 22, 19, 13, 26, 13, 14,
                 22, 24, 21, 22, 26, 21, 23, 24, 27, 41, 31, 27, 35, 26, 28, 36,
                 39, 21, 17, 22, 17, 19, 15, 34, 10, 15, 22, 18, 15, 20, 15, 22,
                 19, 16, 30, 27, 29, 23, 20, 16, 21, 21, 25, 16, 18, 15, 18, 14,
                 10, 15, 8, 15, 6, 11, 8, 7, 18, 16, 13, 12, 13, 20, 15, 16, 12,
                 18, 15, 16, 13, 15, 16, 11, 11)

plot(earthq.data)

model <- depmix(response = earthq.data ~ 1,
                nstates=2,
                family=poisson(),
                data=data.frame(earthq.data = earthq.data))
fitted_model <- fit(model)
fitted_model
# AIC:  693.7574 
# BIC:  707.1216 
summary(fitted_model)

model <- depmix(response = earthq.data ~ 1,
                nstates=3,
                family=poisson(),
                data=data.frame(earthq.data = earthq.data))
fitted_model <- fit(model)
fitted_model
# AIC:  679.055 
# BIC:  708.4561 
summary(fitted_model)

model <- depmix(response = earthq.data ~ 1,
                nstates=4,
                family=poisson(),
                data=data.frame(earthq.data = earthq.data))
fitted_model <- fit(model)
fitted_model
# AIC:  693.1303 
# BIC:  743.9141 
summary(fitted_model)

# nstates should theoretically be 3
posterior(fitted_model)
