using DifferentialEquations
using Plots


# BM model for beta(t): d beta(t) = beta(t) d W(t)
# Function to compute rates of change 
spsir_bm = function(du, u, p, t)
    S = u[1]
    I = u[2]
    R = u[3]
    beta = max(0., u[4])
    
    N = S + I + R
    gamma = p[1] 
    
    du[1] = -beta * S * I / N
    du[2] = beta * S * I / N - gamma * I
    du[3] = gamma * I
    du[4] = 0.
end

# Function to add noise
sigma_spsir = function(du, u, p, t )
    sigma = p[2]
    
    du[1] = 0.
    du[2] = 0.
    du[3] = 0.
    du[4] = sigma 
end

# BM for logbeta(t), with drift 
spsir_logbm_drift = function(du, u, p, t)
    S = u[1]
    I = u[2]
    R = u[3]
    beta = exp( u[4] )
    
    N = S + I + R
    gamma = p[1] 
    alpha = p[3]
    
    du[1] = -beta * S * I / N
    du[2] = beta * S * I / N - gamma * I
    du[3] = gamma * I
    du[4] = -alpha * I
end

# set random seed
srand( 1111 )

## Simulation of BM model 
# starting conditions
u0 = [50.;1.0;0.0;2.0]
tspan = (0.0,10.0)
# parameters gamma sigma 
p = [1.; 1.]
spsir_bm_prob = SDEProblem(spsir_bm, sigma_spsir, u0, tspan, p)
spsir_bm_sol = solve(spsir_bm_prob)

## Simulation of log-BM with drift model 

# starting condtions 
u0 = [50.;1.0;0.0;log(3.) ]
tspan = (0.0,10.0)
# parameters gamma sigma alpha
p = [1.; 1.; 0.1]
spsir_logbm_drift_prob = SDEProblem(spsir_logbm_drift, sigma_spsir, u0, tspan, p)
spsir_logbm_drift_sol = solve(spsir_logbm_drift_prob)

## Plotting for BM model 
# Plot evolution of number infected
plot( spsir_bm_sol , vars = 2 )

# Plot evolution of transmission rate 
plot( spsir_bm_sol , vars = 4 )

## Plotting for BM with drift model 
# Plot evolution of number infected
plot( spsir_logbm_drift_sol , vars = 2 )

# Plot evolution of transmission rate 
plot( spsir_logbm_drift_sol , vars = 4 )
