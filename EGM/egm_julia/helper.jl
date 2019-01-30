struct params
    β::Float64;  #discount factor
    r::Float64;     #interest rate
    σ::Number;      #CRRA
    N::Int64;     #number of gridpoints
end

################################################################################

"""
    Input:  par:        model parameters
            y:          income states
            Π:          income states transition matrix
            grd:        asset grid
            u1:         u'(⋅)
            u1inv:      [u'(⋅)]-¹
            tol:        convergence criterion
            maxiter:    maximum number of iterations
    Output: Cpol_new:   consumption policy function after convergence
            diff:       distance in last update of policy function
Computes the policy function for a simple consumption-savings problem.
"""
function egm(par::params,y::Array{Float64,1},Π::Array{Float64,2},grd::Array{Float64,1},u1,u1inv,tol=1e-6,maxiter=5_000)
    I = length(y);  #number of income states
    N = par.N;      #number of grid points
    β = par.β;
    r = par.r;
    #guesses for policy functions
    Cpol_old = repeat(grd,1,I);   #consumption
    #new policy functions
    Cpol_new = zeros(N,I);
    global ci = zeros(N);  #consumption
    global ai =  zeros(N); #assets

    dif = Inf;
    iter = 0;
    println("\n");
    println("Starting EGM algorithm...");
    while dif>tol && iter<maxiter
        for i=1:I
            ci = (u1inv(β*(1+r)*u1(Cpol_old)*Π[i:i,:]'))[:]; #consumption (c)
            ai = ((ci .+ agrid .- y[i])/(1+r))[:]; #assets today (a)

            #Interpolation
            for n=1:N
                if grd[n]<=ai[1]    #borrowing constraint binds
                    Cpol_new[n,i] = grd[n]*(1+r) + y[i] - grd[1];
                else                #borrowing constraint doesn't bind
                    k= findlast(grd[n].>ai);  #first new asset smaller than grd[n]
                    if k<N-1
                        w = (ai[k+1]-grd[n])/(ai[k+1]-ai[k]); #weight on 1st observation
                        Cpol_new[n,i] = w*ci[k]+(1-w)*ci[k+1];
                    else
                        Cpol_new[n,i] = ci[k];
                    end
                end
            end
        end
        dif = maximum(abs.(Cpol_new.-Cpol_old)[:]);
        if iter%10==0
            println("Iteration #"*string(iter)*", dif = "*string(dif));
        end
        iter+=1;
        Cpol_old = copy(Cpol_new); #updating policy function
    end
    println("Ending EGM.");
    println("\n");
    return Cpol_new, dif;
end

################################################################################

#Plotting function

function plot_fig(agrid,cpol,filename)
    figure()
    for i=1:size(cpol,2)
        plot(agrid,cpol[:,i]);
    end
    xlim([0,30]);
    #ylim([0,2.5]);
    title("Consumption policy functions",fontsize=14);
    xlabel("Savings",fontsize=13);
    ylabel("Consumption",fontsize=13);
    savefig(filename)
end
