#### R Code

## Load packages
library(zoo)
library(tseries)
library(export)
library(MASS)
library(latex2exp)


## Get S&P500 data from Yahoo! Finance
SnP = get.hist.quote("^GSPC", start='2000-01-01',
                     end="2020-12-31", quote="AdjClose")

## Calculate Net Returns
SnP.Ret = diff(SnP)/as.numeric(SnP[-length(SnP)]) # Net returns

## Initial Plot of Index and Net returns
plot(SnP, 
     ylab='S&P 500 Index',
     xlab='Date')

plot(SnP.Ret,
     ylab='S&P 500 Net Returns',
     xlab='Date')

qqnorm(SnP.Ret, datax=TRUE, main='');
qqline(SnP.Ret,datax=TRUE, main='') # Normal qqplot

## Unnormalized density
# takes as input parameters, data and temperature returns value of 
# log of un-normalized density
log_g = function(nu, Y_observed, mu, sigma, temp) {
  
  if (nu<=0 || sigma <= 0){
    return(-Inf) # Reject if v or w <= 0
  }
  
  else {
    n=length(SnP.Ret)
    
    val = (1/temp) * (n*log(gamma((nu+1)/2)/(sigma*sqrt(nu*pi)*gamma(nu/2))) - 
                        ((nu+1)/2) * sum(log(((nu + 
                        ((Y_observed-mu)/sigma)^2))/nu)) + log(nu) -
                        nu/10 - (nu+1)/2 * log(1+sigma^2/2) - mu^2/200)
    # As calculated in previous sections
    
    return(val)
  }
  
}

## ---- Metropolis Algorithm and Simulated Annealing Parameters -----
M = 2*10^5  # run length
temp = 100  # initial value, left constant if no annealing
finaltemp = 0.002 # Final temperature
tempfactor = (temp/finaltemp)^(-1/M)  # for exponential cooling
tempdiff = (temp-finaltemp)/M    # for linearcooling
sigma = 1/3  # proposal scaling
do_tempering=FALSE # Flag for tempering
exp_cool_flag=FALSE # Flag for exponential cooling
seed_val=4 # seed for reproducability

## Function that takes as input above parameters and 
# returns chains of our parameters
compwise_mc = function(M, temp, finaltemp, tempfactor, tempdiff,
                       sigma, do_tempering, exp_cool_flag, seed_val){
  
  set.seed(seed_val) # set seed
  
  # for keeping track of values
  nulist = mulist = siglist = templist = hlist = rep(0,M); numaccept = 0;
  nustart = runif(1, 2, 3) # Overdispersed starting distribution
  mustart = runif(1, 0, 0.002) # Overdispersed starting distribution
  sigstart= runif(1, 0.006, 0.01) # Overdispersed starting distribution
  
  X=c(nustart, mustart, sigstart)
  numxaccept = numtempaccept = temptries = 0;

  # Create Chain, iterate M times
  for (i in 1:M){
    
    coord = sample(1:3, 1) # uniform on {1,3}
    Y = X
    
    if(do_tempering==FALSE){
      temp=1 # No tempering
    }
    
    else{
      
      if(exp_cool_flag==TRUE){
        temp = tempfactor * temp # exponential cooling
      }
      
      else{
        temp = temp - tempdiff;   # linear cooling
      }
      
    }
    
    if (coord==1){
      Y[coord] = X[coord] + sigma * rnorm(1)  # propose in direction "coord"
    }
    else if (coord==2){
      Y[coord] = X[coord] + sigma/800 * rnorm(1)  # propose in direction "coord"
    }
    else{
      Y[coord] = X[coord] + sigma/500 * rnorm(1)  # propose in direction "coord"
      }

    # for accept/reject
    U = log(runif(1))
    alpha = log_g(nu=Y[1], Y_observed=SnP.Ret, mu=Y[2], 
                  sigma = Y[3], temp=temp) - 
    log_g(nu=X[1], Y_observed=SnP.Ret, mu=X[2], sigma = X[3], temp=temp)  

    if (U < alpha){
      X = Y  # accept proposal
      numaccept = numaccept + 1;
      }
      
    nulist[i] = X[1];
    mulist[i] = X[2];
    siglist[i] = X[3];
    templist[i] = temp;
  } 
  cat('Acceptance Rate :', numaccept/M)
  return(list(nulist, mulist, siglist, templist))
}

compwise_metro = compwise_mc(M, temp, finaltemp, tempfactor, tempdiff,
                             sigma, do_tempering, exp_cool_flag, seed_val)

options(scipen=5)

## Degrees of freedom Parameter trace plot
plot(compwise_metro[[1]],type='l', 
     ylab=TeX(r'($\nu_n$)'), xlab="Number of sample in chain (n)", main="")
rolling_mean1=cumsum(compwise_metro[[1]])/(1:M)
lines(rolling_mean1, col='red', lwd=2)

## Degrees of mu Parameter trace plot
plot(compwise_metro[[2]],type='l', 
     ylab=TeX(r'($\mu_n$)'), xlab="Number of sample in chain (n)", main="")
rolling_mean1=cumsum(compwise_metro[[2]])/(1:M)
lines(rolling_mean1, col='red', lwd=2)

## Degrees of sigma Parameter trace plot
plot(compwise_metro[[3]],type='l', 
     ylab=TeX(r'($\sigma_n$)'), xlab="Number of sample in chain (n)", main="")
rolling_mean1=cumsum(compwise_metro[[3]])/(1:M)
lines(rolling_mean1, col='red', lwd=2)

## ACF of parameters
B=20000
mu_df = mean(compwise_metro[[1]][(B+1):M])
acf(compwise_metro[[1]][(B+1):M], lag.max=150, main='')
se1 =  sd(compwise_metro[[1]][(B+1):M]) / sqrt(M-B)
cat("iid standard error would be about", se1, "\n")
varfact <- function(h_x, max_lag){ 2 * sum(acf(h_x, plot=FALSE, 
                                               lag.max=max_lag)$acf) - 1}
vf = varfact(compwise_metro[[1]][(B+1):M], 150)
se = se1 * sqrt( vf )
cat("varfact = ", vf, "\n")
cat("true standard error is about", se, "\n")
cat("approximate 95% confidence interval is (", mu_df - 1.96 * se, ",",
    mu_df + 1.96 * se, ")\n\n")

## Location parameter
mu_df = mean(compwise_metro[[2]][(B+1):M])
acf(compwise_metro[[2]][(B+1):M], lag.max=50, main='')
se1 =  sd(compwise_metro[[2]][(B+1):M]) / sqrt(M-B)
cat("iid standard error would be about", se1, "\n")
varfact <- function(h_x, max_lag){ 2 * sum(acf(h_x, plot=FALSE, 
                                               lag.max=max_lag)$acf) - 1}
vf = varfact(compwise_metro[[2]][(B+1):M], 50)
se = se1 * sqrt( vf )
cat("varfact = ", vf, "\n")
cat("true standard error is about", se, "\n")
cat("approximate 95% confidence interval is (", mu_df - 1.96 * se, ",",
    mu_df + 1.96 * se, ")\n\n")

## Scale parameter
mu_df = mean(compwise_metro[[3]][(B+1):M])
acf(compwise_metro[[3]][(B+1):M], lag.max=150, main='')
se1 =  sd(compwise_metro[[3]][(B+1):M]) / sqrt(M-B)
cat("iid standard error would be about", se1, "\n")
varfact <- function(h_x, max_lag){ 2 * sum(acf(h_x, plot=FALSE, 
                                               lag.max=max_lag)$acf) - 1}
vf = varfact(compwise_metro[[3]][(B+1):M], 150)
se = se1 * sqrt( vf )
cat("varfact = ", vf, "\n")
cat("true standard error is about", se, "\n")
cat("approximate 95% confidence interval is (", mu_df - 1.96 * se, ",",
    mu_df + 1.96 * se, ")\n\n")


## Run algorithm 20 times and see if we get similar results
Nruns=20; B=20000
comp_nu_list=comp_mu_list=comp_sig_list=rep(0, Nruns)
for(i in 1:Nruns){
  vals_compwise = compwise_mc(M, temp, finaltemp, tempfactor, tempdiff,
                              sigma, do_tempering, exp_cool_flag, seed_val=20+i)
  comp_nu_list[i] = mean(vals_compwise[[1]][(B+1):M])
  comp_mu_list[i] = mean(vals_compwise[[2]][(B+1):M])
  comp_sig_list[i] = mean(vals_compwise[[3]][(B+1):M])
  
}

comp_mu_fin = mean(comp_mu_list)
se_fin = sd(comp_mu_list)/sqrt(length(comp_mu_list))
cat("----Sigma= 0.08----\n" ,
    "Final Monte Carlo Estimate for df:", comp_mu_fin, 
    "\nSE : ", se_fin ,
    "\n95% C.I.:  (",
    comp_mu_fin-1.96*se_fin, ",", comp_mu_fin+1.96*se_fin, ")\n",
    sep='')

comp_nu_fin = mean(comp_nu_list)
se_fin = sd(comp_nu_list)/sqrt(length(comp_nu_list))
cat("----Sigma= 0.08----\n" ,
    "Final Monte Carlo Estimate for df:", comp_nu_fin, 
    "\nSE : ", se_fin ,
    "\n95% C.I.:  (",
    comp_nu_fin-1.96*se_fin, ",", comp_nu_fin+1.96*se_fin, ")\n",
    sep='')

comp_sig_fin = mean(comp_sig_list)
se_fin = sd(comp_sig_list)/sqrt(length(comp_sig_list))
cat("----Sigma= 0.08----\n" ,
    "Final Monte Carlo Estimate for df:", comp_sig_fin, 
    "\nSE : ", se_fin ,
    "\n95% C.I.:  (",
    comp_sig_fin-1.96*se_fin, ",", comp_sig_fin+1.96*se_fin, ")\n",
    sep='')

#### Simulated Annealing
## Exponential Cooling
exp_cooling = compwise_mc(M, temp, finaltemp, tempfactor, tempdiff,
                          sigma, do_tempering = TRUE,
                          exp_cool_flag=TRUE, seed_val)

## Degrees of freedom Parameter
plot(exp_cooling[[1]],type='l', 
     ylab=TeX(r'($\nu_n$)'), xlab="Number of sample in chain (n)", main="")

## Degrees of mu Parameter
plot(exp_cooling[[2]],type='l', 
     ylab=TeX(r'($\mu_n$)'), xlab="Number of sample in chain (n)", main="")

## Degrees of sigma Parameter
plot(exp_cooling[[3]],type='l', 
     ylab=TeX(r'($\sigma_n$)'), xlab="Number of sample in chain (n)", main="")

cat('Estimate for df: ', exp_cooling[[1]][M],
    '\n Estimate for mu: ', exp_cooling[[2]][M],
    '\n Estimate for sigma: ', exp_cooling[[3]][M])

## Run algorithm multiple times and see if we get similar results
Nruns=20
exp_nu_list=exp_mu_list=exp_sig_list=rep(0, Nruns)

for(i in 1:Nruns){
  exp_compwise = compwise_mc(M, temp, finaltemp, tempfactor, tempdiff,
                              sigma, do_tempering = TRUE, 
                             exp_cool_flag=TRUE, seed_val=20+i)
  exp_nu_list[i] = exp_compwise[[1]][M]
  exp_mu_list[i] = exp_compwise[[2]][M]
  exp_sig_list[i] = exp_compwise[[3]][M]
  
}

# Standard Deviation of estimates
sd(exp_nu_list)
sd(exp_mu_list)
sd(exp_sig_list)

##### Linear Cooling
linear_cooling = compwise_mc(M, temp, finaltemp, tempfactor, tempdiff,
                             sigma, do_tempering = TRUE, exp_cool_flag=FALSE,
                             seed_val=seed_val)

## Degrees of freedom Parameter
plot(linear_cooling[[1]],type='l', 
     ylab=TeX(r'($\nu_n$)'), xlab="Number of sample in chain (n)", main="")

## Degrees of mu Parameter
plot(linear_cooling[[2]],type='l', 
     ylab=TeX(r'($\mu_n$)'), xlab="Number of sample in chain (n)", main="")

## Degrees of sigma Parameter
plot(linear_cooling[[3]],type='l', 
     ylab=TeX(r'($\sigma_n$)'), xlab="Number of sample in chain (n)", main="")

cat('Estimate for df: ', linear_cooling[[1]][M],
    '\n Estimate for mu: ', linear_cooling[[2]][M],
    '\n Estimate for sigma: ', linear_cooling[[3]][M])

## Run multiple times to see if we get similar results
Nruns=20
lin_nu_list=lin_mu_list=lin_sig_list=rep(0, Nruns)

for(i in 1:Nruns){
  linear_compwise = compwise_mc(M, temp, finaltemp, tempfactor, tempdiff,
                             sigma, do_tempering = TRUE, exp_cool_flag=FALSE, 
                             seed_val = 20+i)
  lin_nu_list[i] = linear_compwise[[1]][M]
  lin_mu_list[i] = linear_compwise[[2]][M]
  lin_sig_list[i] = linear_compwise[[3]][M]
  
}

sd(lin_nu_list)
sd(lin_mu_list)
sd(lin_sig_list)


#### Value at risk calculations

## Generate simulations using our estimated parameters
return_sims_exp = rt(length(SnP.Ret), df=exp_cooling[[1]][M])*
  exp_cooling[[3]][M] + exp_cooling[[2]][M]

return_sims_lin = rt(length(SnP.Ret), df=linear_cooling[[1]][M])*
  linear_cooling[[3]][M] + linear_cooling[[2]][M]

return_sims_metro = rt(length(SnP.Ret), df=comp_nu_fin)*
  comp_sig_fin + comp_mu_fin

plot(density(SnP.Ret), lwd=1.5, main='', 
     xlab='Net Returns', ylab='Kernel Density Estimate')

lines(density(return_sims_exp), col='red', lwd=1.5)
lines(density(return_sims_lin), col='blue',lwd=1.5)
lines(density(return_sims_metro), col='green', lwd=1.5)
text <- c('Historical data', 'Exponential Cooling Sim', 
          'Linear Cooling Sim', 'Metropolis Sim')

legend('topleft', legend=text,
       col=c('black', 'red', 'blue','green'), lwd=1.5, y.intersp=0.8)


alpha=0.05 # Alpha of VaR/CVaR
## Historical VaR
(VaR.historical = -quantile(SnP.Ret,alpha)) # VaR
(CVaR.historical = -mean(SnP.Ret[which(SnP.Ret<(-VaR.historical))])) # CvaR
(exp_ret.historical = quantile(SnP.Ret,0.5)) # Expected Returns

## VaR using Exponential Cooling
t.quant = qt(alpha, exp_cooling[[1]][M])
(VaR.exp_cool = -exp_cooling[[2]][M] - exp_cooling[[3]][M] * t.quant) # VaR
(CVaR.exp_cool = -exp_cooling[[2]][M] + exp_cooling[[3]][M]/alpha *
    dt(t.quant,exp_cooling[[1]][M]) *
    (exp_cooling[[1]][M]+t.quant^2)/(exp_cooling[[1]][M]-1) ) # CVaR
exp_cooling[[2]][M] # Expected returns

## VaR using Linear Cooling
t.quant = qt(alpha, linear_cooling[[1]][M])
(VaR.exp_cool = -linear_cooling[[2]][M] - linear_cooling[[3]][M] * t.quant) # VaR
(CVaR.exp_cool = -linear_cooling[[2]][M] + 
    linear_cooling[[3]][M]/alpha * dt(t.quant,linear_cooling[[1]][M]) *
    (linear_cooling[[1]][M]+t.quant^2)/(linear_cooling[[1]][M]-1) ) # CVaR
linear_cooling[[2]][M] # Expected returns

## VaR using Metropolis Hastings
t.quant = qt(alpha, comp_nu_fin)
(VaR.exp_cool = -comp_mu_fin - comp_sig_fin * t.quant) # VaR
(CVaR.exp_cool = -comp_mu_fin + comp_sig_fin/alpha * dt(t.quant,comp_nu_fin) *
    (comp_nu_fin+t.quant^2)/(comp_nu_fin-1) ) # CVaR
comp_mu_fin # Expected returns
