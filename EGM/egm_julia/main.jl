## Packages
using PyPlot
include("helper.jl");

## Parameters
β = 0.96; #discount factor
r = 0.02; #interest rate
σ = 5;    #EIS
N = 1000;
par = params(β,r,σ,N);
amin = 1e-6;
amax = 100;

#Income process
y = [0.5;4;10];
Π = [0.55 0.4 0.05;
    0.3 0.5 0.2;
    0.1 0.3 0.6];

#Utility
u(x) = (x.^(1-σ))./(1-σ);   #u(⋅)
u1(x) = x.^(-σ);            #u'(⋅)
u1inv(x) = x.^(-1/σ);        #inv(u'(⋅))
#Grid
agrid = e.^(linspace(log(amin),log(amax),N));

#3.1
#egm
cpol, dif = egm(par,y,Π,agrid,u1,u1inv);
plot_fig(agrid,cpol,"EGMpol.png")
close()
