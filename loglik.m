function nll = loglik(mu,sigma,RT,H)

% Returns the negative log likelihood of our model

N = length(RT);
p = normcdf(RT,mu,sigma); %eqn for psychometric fxn

% (google "loglikelihood of psychometric fxn")
for i=1:N
    ll(i,1) = H(i)*log(1/8 + (.95-1/8)*p(i)) + (1-H(i))*log(7/8 - (.95-1/8)*p(i)); %likelihood at each intensity
end

overall_ll = sum(ll);
nll = -overall_ll;

end