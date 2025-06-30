rm(list=ls())     # Clean memory
cat("\014")       # Clear Console
graphics.off()    # Close graphs



#install.packages("nleqslv")      # Package to solve a system of non-linear equations.
#install.packages("maxLik")      # Package to optimize the objective function.
library(nleqslv)
library(maxLik)



##########################################################
################ Parameters specification ################
##########################################################

n <- 1000      # sample size.

set.seed(42)      # to make sure that we always generate the same data.

mu <- 0.5      # true value of mu.

sigma <- 2      # standard deviation of the error term

n.sim <- 500      # number of simulations



###########################################
################ Functions ################
###########################################

# function for the maximum likelihood estimation

ML <- function(par, y){      # objective function called ML; "par" an argument of the function and it represents the vector of the unknown parameters that have to be estimated.
  n <- length(y)
  gamma <- par[1]
  theta <- par[2]
  loglike <- rep(NA,n)      # vector to store the observation-specific log-likelihood values.
  for (i in 1:n){
    logl_discr <- 0      # initialize the value for the discrete distribution "part" of the log-likelihood.
    logl_cont <- 0      # initialize the value for the continuous distribution "part" of the log-likelihood.
    if (y[i] == 0) {
      logl_discr <- log(pnorm(-gamma,0,1))       # expression for the discrete distribution "part" of the log-likelihood.
    } else {
      logl_cont <-  -1/2*(log(2*pi) - log(theta**2) +(theta*y[i] - gamma)**2)     # expression for the continuos distribution "part" of the log-likelihood.
    }
    loglike[i] <- logl_discr + logl_cont      # individual (that is, for observation i) log-likelihood.
  }
  return(loglike)      # a vector of observation-specific log-likelihood values must be returned for the optimization function.
}



# function for the method of moments estimation

MM_equations <- function(par, y) {
  moment_eq <- rep(NA,2)      # vector to store the solutions for the moment equations
  mu.hat <- par[1]
  sigma.hat <- par[2]
  # Define the variables from the derivations
  n_1 <- length(y[y>0])
  n <- length(y)
  lambda_h <- dnorm((-mu.hat/sigma.hat))/(1-pnorm(-mu.hat/sigma.hat))
  # insert the two methods of moments equations here:
  moment_eq[1] <- (n_1/n) - pnorm((mu.hat/sigma.hat))
  moment_eq[2] <- 1/n_1 * sum(y) - mu.hat - sigma.hat
  moment_eq
}



#########################################################
################ Monte Carlo simulations ################
#########################################################

# vectors to store the estimated parameters
mu.hat.ML <- rep(NA,n.sim)
sigma.hat.ML <- rep(NA,n.sim)
mu.hat.MM <- rep(NA,n.sim)
sigma.hat.MM <- rep(NA,n.sim)


# Start the simulations

for (s in 1:n.sim){
  
  y_star <- mu + rnorm(n,0,sigma)      # generate the latent variable.
  
  y <- y_star*(y_star > 0)      # observed variable.
  
  init.val <- c(0.0005,0.0002) 
  
  
  # Maximum likelihood estimation:
  
  # Let's optimize using the Newthon-Raphson algorithm:
  optNR <- maxNR(ML, start=init.val, y=y, iterlim=20)      # Maximize the objective function with the Newthon-Raphson method.
  par.hat <- optNR$estimate      # the ML estimates for gamma and theta.
  mu.hat.ML[s] <-   par.hat[1]/par.hat[2]    # the ML estimate for mu.
  sigma.hat.ML[s] <- 1/par.hat[2]       # the ML estimate for sigma.
  
  
  # Method of moments estimation:
  
  result <- nleqslv(init.val, MM_equations, y=y)
  mu.hat.MM[s] <- result$x[1] 
  sigma.hat.MM[s] <- result$x[2] 
  
}


# Compute the mean squared errors for each estimator
MSE.mu.hat.ML= mean((mu.hat.ML-mu)^2)
MSE.sigma.hat.ML= mean((sigma.hat.ML-sigma)^2)
MSE.mu.hat.MM= mean((mu.hat.MM-mu)^2)
MSE.sigma.hat.MM= mean((sigma.hat.MM-sigma)^2)
# outputs for this instance with n=1000:
# > MSE.mu.hat.ML
#[1] 0.004834958
#> MSE.sigma.hat.ML
#[1] 0.003971676
#> MSE.mu.hat.MM
#[1] 0.08212364
#> MSE.sigma.hat.MM
#[1] 1.306659
#Clearly MLE is better, which should make sense, since
# MLE is efficient, so this should give the better estimates (at least quicker)


# print results to screen
MSE.mu.hat.ML
MSE.sigma.hat.ML
MSE.mu.hat.MM
MSE.sigma.hat.MM







import pandas as pd

# ====== PARAMETRY ======

file_path = "twoj_plik.xlsx"

# Grupy walut
G4 = ["EUR", "USD", "GBP", "JPY"]
Other_G10 = ["AUD", "CAD", "NZD", "NOK", "SEK", "CHF"]
Other_non_G10 = ["BRL", "CNY", "DKK", "HKD", "KRW", "MXN", "RUB", "SGD", "TRY", "ZAR"]

# ====== WCZYTYWANIE SHEETÓW ======

df_pln = pd.read_excel(file_path, sheet_name="PLN")
df_usd = pd.read_excel(file_path, sheet_name="USD")
df_eur = pd.read_excel(file_path, sheet_name="EUR")
df_other = pd.read_excel(file_path, sheet_name="Others")

# ====== AGREGOWANIE SUM PO WALUTACH ======

def sum_columns(df, exclude_cols=["others", "residuals"]):
    cols = [col for col in df.columns if col not in exclude_cols]
    sums = df[cols].sum()
    return sums

# Sumy z jawnych kolumn
sums_pln = sum_columns(df_pln)
sums_usd = sum_columns(df_usd)
sums_eur = sum_columns(df_eur)

# Sumy z ostatniego sheeta (rozbite others + residuals)
sums_other = df_other.sum()

# ====== ŁĄCZENIE SUM PO WALUTACH ======

total_sums = pd.concat([sums_pln, sums_usd, sums_eur, sums_other])
total_sums = total_sums.groupby(level=0).sum()

# ====== SUMOWANIE PER GRUPA ======

group_sums = {"PLN": 0, "G4": 0, "Other_G10": 0, "Other_non_G10": 0, "Other": 0}

for currency, value in total_sums.items():
    if currency == "PLN":
        group_sums["PLN"] += value
    elif currency in G4:
        group_sums["G4"] += value
    elif currency in Other_G10:
        group_sums["Other_G10"] += value
    elif currency in Other_non_G10:
        group_sums["Other_non_G10"] += value
    else:
        group_sums["Other"] += value

# ====== TOTAL PLN + TOTAL ALL ======

total_pln = group_sums["PLN"]
total_all = sum(group_sums.values())

# ====== WYŚWIETLENIE PODSUMOWANIA ======

print("✅ Podsumowanie per grupa:")
for group, val in group_sums.items():
    print(f"{group}: {val}")

print(f"\n🔹 TOTAL PLN: {total_pln}")
print(f"🔹 TOTAL ALL CURRENCIES: {total_all}")

print("\n✅ Podsumowanie per waluta:")
print(total_sums)

# ====== ZAPIS DO EXCELA ======

summary_df = pd.DataFrame(list(group_sums.items()), columns=["Group", "Total"])
total_sums_df = total_sums.reset_index()
total_sums_df.columns = ["Currency", "Total"]

totals_df = pd.DataFrame({
    "Metric": ["Total PLN", "Total All Currencies"],
    "Value": [total_pln, total_all]
})

with pd.ExcelWriter("podsumowanie.xlsx") as writer:
    summary_df.to_excel(writer, sheet_name="Group_Summary", index=False)
    total_sums_df.to_excel(writer, sheet_name="Currency_Summary", index=False)
    totals_df.to_excel(writer, sheet_name="Totals", index=False)

print("\n✅ Plik 'podsumowanie.xlsx' został utworzony!")
