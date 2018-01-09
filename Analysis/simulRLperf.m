function perfRL = simulRLperf(blck,n,alpha,beta)

if nargin==1
    n = blck.ntrl;
    alpha = .3;
    beta  = 5;
elseif nargin == 2
    alpha = .3;
    beta  = 5;    
end

deltaQ = deltaQlearn(alpha,blck.outcome(1:n));

actions = ((rand(n,1))<(1./(1+exp(beta*deltaQ))))*1.;
actions(actions == 0) = 2;

perfRL = sum(actions' == blck.outcome(1:n))/n;

end