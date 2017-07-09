## example:

using Distributions, Gadfly

#setting parameters
ϕ₁ = 0.5;
ϕ₂ = -0.2;
ϕ₃ = 0;
ϕ₄ = 0.5;
σ = 0.2;
T = 100;
n_tries = 500000;

## drawing randomness
d = Normal();

x0 = [1;1;1;1];


G = [1 0 0 0];
A = [ϕ₁ ϕ₂ ϕ₃ ϕ₄;
     1 0 0 0;
     0 1 0 0;
     0 0 1 0];
C = [σ; 0; 0; 0];

YT = zeros(n_tries,1);

for n=1:n_tries;
    w = rand(d,T); #drawng noise
    y = x0; #setting start at x0
    for t=1:T;
        tmp = A*y + C*w[t];
        YT[n] = (G*tmp)[1];
        y = tmp;
    end
end

plot(x=YT,Geom.histogram(bincount=100))
