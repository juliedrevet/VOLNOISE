function [a1,b1,a2,b2,a3,b3] = RLfit_vol(expe)

nblck = length(expe.blck);
nvol = length(expe.blck(1).volperm);
ntrl  = length(expe.blck(1).outcome)/nvol;
outcome = zeros(nvol,ntrl*nblck);
choice  = zeros(nvol,ntrl*nblck);

for iblck = 1:nblck
    for ivol = 1:nvol
        idx_vol = (expe.blck(iblck).volperm(ivol));
        outcome(idx_vol,((iblck-1)*ntrl+1):iblck*ntrl) = expe.blck(iblck).outcome(((ivol-1)*ntrl+1):ivol*ntrl);
        choice(idx_vol,((iblck-1)*ntrl+1):iblck*ntrl) = expe.rslt(iblck).resp(((ivol-1)*ntrl+1):ivol*ntrl);
    end
end

Eq1=@(alpha,beta)(neglog(deltaQlearn(alpha,outcome(1,:),ntrl),choice(1,:),beta));
LLH1=@(X)Eq1(X(1),X(2));
Eq2=@(alpha,beta)(neglog(deltaQlearn(alpha,outcome(2,:),ntrl),choice(2,:),beta));
LLH2=@(X)Eq2(X(1),X(2));
Eq3=@(alpha,beta)(neglog(deltaQlearn(alpha,outcome(3,:),ntrl),choice(3,:),beta));
LLH3=@(X)Eq3(X(1),X(2));

alphabeta1 = fminsearch(LLH1,[0.01,5]);
alphabeta2 = fminsearch(LLH2,[0.01,5]);
alphabeta3 = fminsearch(LLH3,[0.01,5]);

a1 = alphabeta1(1); b1 = alphabeta1(2);
a2 = alphabeta2(1); b2 = alphabeta2(2);
a3 = alphabeta3(1); b3 = alphabeta3(2);
end