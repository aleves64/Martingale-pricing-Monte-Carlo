S0 <- 50 # Spot price
sigma <- 0.3 # Volatility
T <- 1 # Time to maturity
dt <- 1/252 # Time step
m <- 10**5 # Number of simulations
r <- 0.05 # Annual risk free interest rate
K <- 50 # Strike price

# Random variates
eps <- rnorm(m)
# Stock price at time T under martingale measure with exp(r*t) as the numeraire
S_T_Q <- S0 * exp((r-0.5*sigma**2)*T + sigma*eps*sqrt(T))
# Stock price at time T under martingale measure with stock price as the numeraire
S_T_Qs <- S0 * exp((r + sigma**2 - 1/2 * sigma**2)*T + sigma*eps*sqrt(T))

# Price under Q
C0_MC_Q <- exp(-r*T)*mean(pmax(S_T_Q - K, 0))
# Price under Qs
C0_MC_Qs <- S0*mean(pmax((S_T_Qs - K)/S_T_Qs, 0))
# Price using both
C0_MC_Both <- S0*mean(S_T_Qs > K) - K*exp(-r*T)*mean(S_T_Q > K)
# Price with Black-Scholes
d1 <- (log(S0/K)+(r+0.5*sigma**2)*T)/(sigma*sqrt(T))
d2 <- d1 - sigma*sqrt(T)
C0_BS <- S0*pnorm(d1) - K*exp(-r*T)*pnorm(d2)

results <- data.frame(Price_under_Q=C0_MC_Q,
                      Price_under_Qs=C0_MC_Qs,
                      Price_using_both=C0_MC_Both,
                      Black_Scholes_price=C0_BS)
results
