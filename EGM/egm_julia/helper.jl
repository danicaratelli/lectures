"""
    Input:  ρ:          autoregressive coefficient of log-income
            sig_y:      log-income standard deviation
            m:          number by which to approximate with discrete distribution
            multi:      number of stds for discrete support (default = 3)
    Output: ygrid:    residual standard deviation
            Pi:         transition matrix (Markov)
            p:         stationary distribution
Approximation of continuous log-income process with discrete process.
Distributions is a required package
"""
function tauchen_income(ρ::Float64,sig_y::Float64,m::Int,multi::Int=3)
    sig_eps = sig_y*sqrt(1-ρ^2);
    dNorm = Normal(0,sig_eps);  #distribution of residuals
    ygrid = linspace(-multi*sig_y,multi*sig_y,m);
    d = diff(ygrid)[1]/2;

    #computing transition matrix
    Pi = zeros(m,m);
    for j=2:m-1
        eps1 = ygrid[1] .+d .-ρ*ygrid; #residuals going to j location
        Pi[:,1] = cdf(dNorm,eps1);
        epsj = ygrid[j] .-ρ*ygrid; #residuals going to j location
        Pi[:,j] = cdf(dNorm,epsj.+d) - cdf(dNorm,epsj.-d);
    end
    epsm = ygrid[m] .-d .-ρ*ygrid; #residuals going to j location
    Pi[:,m] = 1-cdf(dNorm,epsm);

    #computing stationary distribution
    tol = Inf;
    p = (1/m)*ones(1,m);
    while tol>1e-5
        p_new = p*Pi;
        tol = maximum(abs.(p_new-p));
        p = p_new
    end

    #normalizing income to 1
    y = exp.(ygrid)/(sum(p*exp.(ygrid)));
    return y, Pi, p;
end

################################################################################

"""
    Input:  c:      consumption
            agrid:  asset grid
            ygrid:  income grid
            r:      interest rate
    Output: a:      new asset grid
            Indxs:  indexes where budget constraint does not bind
            W:      weights from interpolation
Finds asset grid consistent with next period assets and current consumption for
each (current) income state.
"""
function interp(c::Array{Float64,2},agrid::Array{Float64,1},y::Array{Float64,1},r::Float64)
    m,N = size(c);

    a = zeros(m,N);     #current assets
    W = zeros(m,N);     #weights
    Indxs = zeros(m,N); #first index where LHS<=RHS
    for i=1:m
        LHS = c[i,:] + agrid;           #c + a
        RHS = (1+r)*agrid + y[i];   #(1+r)*a + y (cash on hand)
        for j=1:N
            rhsj = RHS[j];
            j_index = findlast(LHS.<rhsj);
            if j_index.==nothing #budget constraint does not hold, i.e. LHS > RHS
                j_index = 1;
            elseif j_index>1     #budget constraint holds, i.e. LHS <= RHS
                j_index-=1;
            end
            #linear interpolation
            w = (LHS[j_index + 1] - rhsj) / (LHS[j_index + 1] - LHS[j_index]);  #weight for linear interpolation
            a[i,j] = w*agrid[j_index] + (1 - w)*agrid[j_index+1];   #interpolation of asset point

            W[i,j] = w;
            Indxs[j] = j_index;
        end
    end
    return a, Indxs, W
end
