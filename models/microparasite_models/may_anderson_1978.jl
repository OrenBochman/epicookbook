using DifferentialEquations
using Plots


function macrop_ode(dY,Y,p,t)
 dY[1] = (p[1]-p[2])*Y[1] - p[3]*Y[2]
 dY[2] = p[4]*Y[1]*Y[3] - (p[5]+p[3]+p[2])*Y[2] - (p[3]*((Y[2]^2)/Y[1])*((p[6]+1)/p[6])) 
 dY[3] = p[7]*Y[2] - (p[8]*Y[3]) - (p[4]*Y[1]*Y[3])
end

par=[1.4,1.05,0.0003,0.01,0.5,0.1,10.0,10.0]
init=[100.0,10.0,10.0]
tspan=(0.0,100.0)

macro_par = ODEProblem(macrop_ode,init,tspan,par)

sol=solve(macro_par);

# Plot
plot(sol,xlabel="Time",yscale=:log10)
