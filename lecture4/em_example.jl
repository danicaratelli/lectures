## Here we have an example of the application of the EM-algorithm, following the example given in What is the expectation maximization algorithm? by Do and Batzoglou

#=
the example is that of two (possibly biased) coins of which we want to estimate the bias, more precisely, the number of heads over all the flips of coin 1 and 2, call this ratio θ₁=(#of heads using 1)/(#flips using 1) and similarly for 2. Then θ = (θ₁,θ₂). We assume that, we do not know which coin is being flipped at each iteration
=#

#preparation
using Distributions
num_iter1 = 500;
num_iter2 = 700;
num_draws = 10;
theta_0 = [0.6,0.5]; #initial guess to the biases
treshold = 0.0001; #treshold to stop iteration

#constructin sample data
theta_1 = 0.35; #actual probability imposed
theta_2 = 0.60; #actual probability imposed



#constructing real data
num_iter = num_iter1 + num_iter2;
d1 = Binomial(1,theta_1);
draws_1 = rand(d1,num_draws*num_iter1);
draws_1 = reshape(draws_1,num_draws,num_iter1);
d2 = Binomial(1,theta_2);
draws_2 = rand(d2,num_draws*num_iter2);
draws_2 = reshape(draws_2,num_draws,num_iter2);
draws = [draws_1 draws_2];


heads = sum(draws,1);
#heads = [5 9 8 4 7] 
#initial conditions
θ = zeros(1,2);
θ[1,:] = theta_0;


#iterative process
j=1;
tr=1;
loglik_old = -Inf;

function lnlik(ps,theta,H,draws);
    tot=0;
    j=1;
    while j<=size(ps,1);
        tot = tot + ps[j,1]*log(theta[1]^H[j]*(1-theta[1])^(draws-H[j]))+
        ps[j,2]*log(theta[2]^H[j]*(1-theta[2])^(draws-H[j]));
        j=j+1;
    end
    return tot;
end

    
while tr>treshold
    #E-step
    #computing likelihoods for each iteration, for each die and determining the probabilities the draws were from each die
    ps = zeros(num_iter,2);
    exps_heads = zeros(num_iter,2);
    exps_tails = zeros(num_iter,2);
    for it = 1:num_iter;
        lik1 = (θ[j,1]^heads[it])*((1-θ[j,1])^(num_draws-heads[it]));
        lik2 = (θ[j,2]^heads[it])*((1-θ[j,2])^(num_draws-heads[it]));
        ps[it,1] = (lik1/(lik1+lik2)); #probability the draws came from coin 1
        ps[it,2] = (lik2/(lik1+lik2)); #probability the draws came from coin 2
        exps_heads[it,1] = ps[it,1]*heads[it]; #expected number of heads in X1 coming from coin 1
        exps_heads[it,2] = ps[it,2]*heads[it]; #expected number of heads in X1 coming from coin 2
        exps_tails[it,1] = ps[it,1]*(num_draws-heads[it]); #expected number oftails in X1 coming from coin 1
        exps_tails[it,2] = ps[it,2]*(num_draws-heads[it]); #expected number of tails in X1 coming from coin 2
    end

    #M-step
    j = j+1;
    theta_new = zeros(1,2);
    theta_new[1] = sum(exps_heads[:,1])/(sum(exps_heads[:,1]+exps_tails[:,1]));
    theta_new[2] = sum(exps_heads[:,2])/sum((exps_heads[:,2]+exps_tails[:,2]));
    θ = [θ;theta_new];
    
    #computing the loglikelihoods
    loglik = lnlik(ps,theta_new,heads,num_draws);
    tr = loglik-loglik_old;
    println("log-likelihood improvement is: "*string(tr));
    loglik_old = loglik;
end


println("your EM-estimated θ is"*string(θ[end,:]))
