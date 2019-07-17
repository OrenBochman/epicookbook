using DifferentialEquations
using IterableTables, DataFrames
using NamedTuples
using Gadfly


function F(du,u,p,t)
    S_H, E_H, I_H, R_H, S_V, E_V, I_V = u
    
    # host dynamics
    host_infection = (p.β*S_H*I_V)/p.N_H
    host_mortality = p.μ_H .* u[1:4] # include S_H, so easier to remove mortality
    host_births = sum(host_mortality)
    host_progression = p.σ_H*E_H
    recovery = p.λ*I_H
    
    du[1] = -host_infection + host_births
    du[2] = host_infection - host_progression
    du[3] = host_progression - recovery
    du[4] = recovery
    du[1:4] -= host_mortality
    
    # vector dynamics
    vec_infection = (p.β*S_V*I_H)/p.N_H
    vec_mortality = p.μ_V .* u[5:7] # include S_V, so easier to remove mortality
    vec_births = sum(vec_mortality)
    vec_progression = p.σ_V*E_V
    
    du[5] = -vec_infection + vec_births
    du[6] = vec_infection - vec_progression
    du[7] = vec_progression
    du[5:7] -= vec_mortality
    
end

# nb: in >= Julia v0.7, can eliminate this import
#  and the @NT syntax
u0 = [
    S_H=100.0,   E_H=0.0, I_H=1.0, R_H=0.0,
    S_V=10000.0, E_V=0.0, I_V=0.0
]
p = @NT(
  μ_H=1/365, μ_V=1/30, σ_H=1/3, σ_V=1/7, λ=1/14,
  β=0.05, N_H = sum(u0[1:4])
)
tspan = (0.0, 365.0)
prob = ODEProblem(F, u0, tspan, p)
sol = @time solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8,saveat=linspace(0,365,365*10+1))

df = DataFrame(sol)
rename!(df, Dict(:timestamp => :t,
  :value1 => :S_H, :value2 => :E_H, :value3 => :I_H, :value4 => :R_H,
  :value5 => :S_V, :value6 => :E_V, :value7 => :I_V
))
mlt = melt(df,:t) # convert results into long format for plotting
mlt[:host] = contains.(string.(mlt[:variable]),"H"); # tag which entries are host vs vector
df

# Plot
fig1a = plot(mlt[mlt[:host] .== true,:], x=:t, y=:value, color=:variable, Geom.line)
fig1b = plot(mlt[mlt[:host] .!= true,:], x=:t, y=:value, color=:variable, Geom.line)
vstack(fig1a,fig1b)
