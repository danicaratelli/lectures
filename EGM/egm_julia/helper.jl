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
    while tol>1e-10
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
            coh:    cash on hand
            r:      interest rate
    Output: a:      new asset grid
            Indxs:  indexes where budget constraint does not bind
            W:      weights from interpolation
Finds asset grid consistent with next period assets and current consumption for
each (current) income state.
"""
function interp(c::Array{Float64,2},agrid::Array{Float64,1},coh::Array{Float64,2},r::Float64)
    m,N = size(c);

    a = zeros(m,N);     #current assets
    W = zeros(m,N);     #weights
    Indxs = zeros(m,N); #first index where LHS<=RHS
    for i=1:m
        LHS = c[i,:] + agrid;           #c + a
        RHS = coh[i,:];   #(1+r)*a + y (cash on hand)
        jindex = 1;
        for j=1:N
            rhsj = RHS[j];
            #jindex first index where budget constraint holds
            while jindex < N - 1
                if LHS[jindex] >= rhsj
                    break
                end
                jindex += 1
            end
            if jindex!==1
                jindex_loc = jindex-1;
            else
                jindex_loc = 1;
            end
            #linear interpolation
            w = (LHS[jindex_loc + 1] - rhsj) / (LHS[jindex_loc + 1] - LHS[jindex_loc]);  #weight for linear interpolation
            a[i,j] = w*agrid[jindex_loc] + (1 - w)*agrid[jindex_loc+1];   #interpolation of asset point

            W[i,j] = w;
            Indxs[j] = jindex_loc;
        end
    end
    return a, Indxs, W
end


################################################################################

"""
    Input:  agrid:  asset grid
            coh:    cash on hand
            c_pol:  consumption policy function
            r:      interest rate
            β:      discount factor
            Pi:     transition matrix
            σ:      constant of risk aversion
    Output: a:      new asset grid
Finds asset grid consistent with next period assets and current consumption for
each (current) income state.
"""
function back_iterate(agrid::Array{Float64,1},coh::Array{Float64,2},c_pol::Array{Float64,2},r::Float64,β::Float64,Pi::Array{Float64,2},σ::Float64)
    c = (β*(1+r)*Pi*((c_pol).^(-σ))).^(-1/σ);
    a, Indxs, W = interp(c,agrid,coh,r);
    a[a.<agrid[1]] = agrid[1];
    return a;
end
