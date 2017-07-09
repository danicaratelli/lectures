using Plots, MAT, Distributions

#importing data
include(joinpath(pwd(),"data","initialize.jl"));

###### running ######
#options
num_pc = 4;


#function setting
function norm_data(X::Array{Float64,2})
    sz = size(X);
    return (X-repmat(mean(X,1),sz[1],1))./repmat(std(X,1),sz[1],1);
end

function pca(X::Array{Float64,2},n::Integer)
    Ω = cov(X); #variance covariance matrix
    E = eig(Ω); eval = E[1]; evec = E[2];
    if n<=size(X,2)
        FV = evec[:,end-n+1:end];
    else
        FV = evec;
    end
      return  PCA_X = FV'*X';
end

function invv(x);
    if ~isempty(size(x));
        return inv(x);
    else
        return 1/x;
    end
end


#assume a normal prior

#function filtering_dist()
#end
    
#organizing data
i=0;
eval(parse("data_tmp = data"*string(i)));
data_tmp = data_tmp.Data;
to_keep = find(sum(isnan(data_tmp),2).==0);
data_tmp = data_tmp[to_keep,:]; #data for which we have all data
sz_new = size(data_tmp); 
data_norm = norm_data(data_tmp);
pca_res = pca(data_norm[:,2:end],num_pc);
β = inv(pca_res*pca_res')*pca_res*data_norm[:,1];
eps = data_norm[:,1] - pca_res'*β;
prior = Normal(0.0, 1);
for i=1:3;
    eval(parse("data_tmp = data"*string(i)));
    data_tmp = data_tmp.Data;
    to_keep = find(sum(isnan(data_tmp),2).==0);
    data_tmp = data_tmp[to_keep,:]; #data for which we have all data
    sz_new = size(data_tmp); 
    data_norm = norm_data(data_tmp);
    pca_res = pca(data_norm[:,2:end],num_pc);

    y = data_norm[end];
    lk = Normal(y,var(eps));

    post_mean = prior.μ+(prior.σ*1*invv(prior.σ + lk.σ)*
end
