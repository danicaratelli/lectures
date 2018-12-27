## Packages
using Distributions, PyPlot
include("helper.jl");

## Parameters
β = 0.92;       #discount factor
σ = 1.0;          #risk aversion parameters
sig_y = 0.92;      #income standard distribution
amax = 10_000;  #maximal assets
r = 0.02;       #interest rate
N = 500;        #number of grid points
m = 11;         #discretization of shock
ρ = 0.9;       #autoregressive coefficient

u(c) = (c.^(1-σ))/(1-σ); #utility function CRRA
uc(c) = c.^(-σ); #marginal utility of consumption

#constructing (next period) asset grid
agrid = linspace(log(0.25),log(amax+0.25),N);
agrid = exp.(agrid) - 0.25;

#discretizing income shock using tauchen.
#Income process is as: log(y_t) = ρlog(y_t-1) + ε
y, Pi, p = tauchen_income(ρ,sig_y,m);

coh = repeat(y,1,N).+ (1+r)*repeat(agrid',m,1) #consumption today c(y,a)

## Section 1 --- 1st iterations
c_pol = coh .+ 0.02*repeat(agrid',m,1);   #guess consumption policy function
c = (β*(1+r)*Pi*uc(c_pol)).^(-1/σ); #consumption today c(y,a') from Euler Equations

#computing assets today (via interpolation)
a, Indxs, W = interp(c,agrid,coh,r);
a[a.<agrid[1]] = agrid[1];  #replacing where assets are too low
#updating consumption policy function (using budget constraint)
c_pol = coh .- a; #budget constraint is c_pol + a' = coh = (1+r)a + y

## Section 2 --- iterate until convergence
epsilon = 1e-6;     #tolerance criterion
tol = Inf;          #current tolerance
max_iter = 1_000;   #maximum number of iterations
it = 1;             #iteration
a_pol_old = copy(a);
a_pol = copy(a_pol_old);
while tol>epsilon && it<max_iter
    if it % 10 == 0
        println("Iteration # "*string(it)*"\n");
        println("Tolerance: "*string(tol));
    end
    a_pol = back_iterate(agrid,coh,c_pol,r,β,Pi,σ);
    c_pol = coh .- a_pol;
    tol = maximum(abs.(a_pol.-a_pol_old));
    a_pol_old = copy(a_pol);
    it+=1;
end
