function [alphabeta, alphabeta_opt] = RLfit(expe)

nblck = length(expe.blck);
ntrl = length(expe.blck(1).outcome);
outcome = zeros(1,ntrl*nblck);
choice  = zeros(1,ntrl*nblck);
correct = zeros(1,ntrl*nblck);

for iblck = 1:nblck
    outcome(((iblck-1)*ntrl+1):iblck*ntrl) = expe.blck(iblck).outcome;
    choice(((iblck-1)*ntrl+1):iblck*ntrl)  = expe.rslt(iblck).resp;
    correct(((iblck-1)*ntrl+1):iblck*ntrl) = expe.blck(iblck).correct;
end

Eq  = @(alpha,beta)(neglog(deltaQlearn(alpha,outcome,ntrl),choice,beta));
LLH = @(X)Eq(X(1),X(2));

Eq_opt  = @(alpha,beta)(neglog(deltaQlearn(alpha,outcome,ntrl),correct,beta));
LLH_opt = @(X)Eq_opt(X(1),X(2));

alphabeta = fminsearch(LLH,[0.01,5]);
alphabeta_opt = fminsearch(LLH_opt,[0.01,5]);
end