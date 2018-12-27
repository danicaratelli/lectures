## Packages
using Distributions, PyPlot
include("helper.jl");

## Parameters
β = 0.92;       #discount factor
σ = 1;          #risk aversion parameters
sig_y = 0.92;      #income standard distribution
amax = 10_000;  #maximal assets
r = 0.02;       #interest rate
N = 500;        #number of grid points
m = 11;         #discretization of shock
ρ = 0.9;       #autoregressive coefficient

u(c) = (c.^(1-σ))/(1-σ); #utility function CRRA
uc(c) = c.^(-σ); #marginal utility of consumption

agrid = linspace(log(0.25),log(amax+0.25),N); #constructing asset grid
agrid = exp.(agrid) - 0.25;

#discretizing income shock using tauchen.
#Income process is as: log(y_t) = ρlog(y_t-1) + ε
y, Pi, p = tauchen_income(ρ,sig_y,m);

coh = repeat(y,1,N).+ (1+r)*repeat(agrid',m,1) #consumption today c(y,a)

#c_pol = coh;
c_pol = coh .+ 0.02*repeat(agrid',m,1);   #guess consumption policy function
c = (β*(1+r)*Pi*uc(c_pol)).^(-1/σ); #consumption today c(y,a') from Euler Equations

#computing assets (via interpolation)
a = interp(c,agrid,y,r);
