# This code accompanies the article XXXX
# It is written for R version 3.6.1, using the package 'data.table'
# It shows how the SMR and life expectancies reported in the article are calculated
# The results vary slightly from the published results because:
#   (a) noise has been added to the OATS data
#   (b) the actual analysis stratifies mortality rates by year as well as age and sex
# Data and results are only provided for men, because the majority of the sample is male and confidentiality was simpler to maintain

library(data.table)
set.seed(17)

#--------------------
# read mortality data
#--------------------

# OATS mortality data. 'fut' = follow-up time (person-years)
# follow-up time has been calculating by age during follow-up (rather than age at baseline), by day
# numbers of deaths have a random value of {rpois(N, 5) - 5} added

m <- fread('https://raw.githubusercontent.com/danlewer/oats/master/oats_mortality_data.csv') 

#--------------------
# read reference data
#--------------------

# this data is the total deaths and aggregate population (i.e. an estimate of person-years) 2001-2018
# reference data is avaiable at:
# Australian Bureau of Statistics. ABS.Stat. 2019. Available from: http://stat.data.abs.gov.au/

r <- fread('https://raw.githubusercontent.com/danlewer/oats/master/ref_nsw_data.csv')

#----------
# join data
#----------

mr <- m[r, on = 'age']

#--------------
# calculate SMR
#--------------

# add age groups
mr[, age_group := findInterval(age, c(18, seq(25, 75, 10)))]
ag <- aggregate(age ~ age_group, mr, range)
ag$lab <- paste0(ag$age[,1], '-', ag$age[,2])
mr[, age_group := factor(age_group, ag$age_group, ag$lab)]

# add expected deaths
mr[, ref_rate := ref_m / ref_p]
mr[, expected := fut * ref_rate]

# summarise by age group
smr <- mr[age %in% 18:74, lapply(.SD, sum), 'age_group'][, -'age']
smr <- rbind(smr, t(c(age_group = NA, colSums(smr[,-1])))) # add total row

# calculate smr
smr[, smr := m / expected]

# monte-carlo confidence intervals
B <- 10000 # number of simulations
mmc <- rpois(B * nrow(smr), smr$m)
mmc <- matrix(mmc, ncol = B) # B simulations
mmc <- mmc / smr$expected # smr values
mmc <- apply(mmc, 1, quantile, prob = c(0.025, 0.5, 0.975))
smr <- cbind(smr, t(mmc))

# note similarity of CI's to classic method (1.96 +/- sqrt(observed) / expected)
# SMR for men = 5.12 (4.96 to 5.29). Value in article is 5.14 (4.97 to 5.30)

#------------------------------------
# mortality at youngest / oldest ages
#------------------------------------

# predict mortality rate based on exponential increase
model <- glm(m ~ age + offset(log(fut)), mr[age %in% 18:74], family = 'poisson')
mr[, pred_rate := predict(model, newdata = data.table(age = 0:99, fut = 100000), type = 'response')]

# 'base' scenario is actual rate from ages 25-64, and maximum of predicted rate from poisson model and general population rate
mr[, oat_rate := m / fut]
mr[, base := ifelse(age %in% 25:64, oat_rate, pred_rate / 100000)]
mr[, base := pmax(ref_rate, base)]

# plot 
plot(1, type = 'n', xlim = c(0, 100), ylim = log(c(10, 100000)), axes = F, xlab = 'Age', ylab = 'Mortality rate / 100000')
lines(mr$age, log(mr$ref_rate * 100000), lty = 3) # general population mortality rate
points(mr$age, log(mr$oat_rate * 100000)) # actual rate in cohort
with(mr[age > 17], lines(age, log(pred_rate))) # modelled rate
with(mr[age > 17], lines(age, log(base * 100000), col = 'red', lty = 3, lwd = 2)) # base case estimate
axis(1, pos = log(10))
axis(2, log(10^(1:5)), 10^(1:5), las = 2, pos = 0)

#--------------------------------
# calculate life expectancy at 18
#--------------------------------

# life table function

life.table <- function(mx, cohort = 100000) {
  n <- length(mx) + 1
  qx <- 2 * mx / (2 + mx)
  qx <- c(qx, 1) # forced method - mortality rate max age + 1 is 100%
  lx <- c(1, cumprod(1 - qx)) * cohort
  dx <- -c(diff(lx), lx[n] * qx[n])
  t <- (lx + c(lx[-1], 0)) / 2
  Tx <- rev(cumsum(rev(t)))
  ex <- Tx / lx
  data.frame(lx = lx, dx = dx, t = t, Tx = Tx, ex = ex)[1:n,]
}

# make life tables for OAT cohort and general population

lt_oat <- life.table(mr[age > 17]$base)
lt_ref <- life.table(mr[age > 17]$ref_rate)

lt_oat$ex[1] # life expectancy for men age 18 = 48.3
lt_ref$ex[1] # general population = 62.2
lt_ref$ex[1] - lt_oat$ex[1] # gap = 13.9 years

# plot survival curves

plot(1, type = 'n', xlim = c(18, 100), ylim = c(0, 100000), xlab = 'Age', ylab = 'Number in cohort surviving')
ages <- 18:100
lines(ages, lt_ref$lx, lty = 3) # general population
lines(ages, lt_oat$lx, lwd = 2, col = 'red') # OAT cohort
