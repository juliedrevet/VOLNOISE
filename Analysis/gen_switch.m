%generate switch sequence
function [y, evol] = gen_switch(n,ntrl,v,csrt,minstd)
% output: switch sequence of ntrl trials for volatility v with contraint csrt
%
% generate n switch sequences and takes the one with smallest zscore
% according to the mean length of episodes also controlling the std of the
% episode lengths

seq  = zeros(n,ntrl);
ml   = zeros(n,1); % mean length before switch
stdl = zeros(n,1); % std length before switch

for m = 1:n
    % generate reward switch sequence
    seq(m,:) = seqswitch(v,csrt,ntrl);
    
    idx = find(seq(m,:));
    
    if isempty(idx)       % nswitch = 0: length is ntrl
        ml(m) = ntrl;
    elseif length(idx)==1 % nswitch = 1: length is the mean of both episodes
        ml(m) = mean([idx, ntrl-idx]);
        stdl(m) = std([idx,ntrl-idx]);
    else                  % nswitch >=2: mean length taking beginning and end of block into account
        ml(m) = mean(diff([1 idx ntrl]));
        stdl(m) = std(diff([1 idx ntrl]));
    end
end

zs = zscore(ml);

idx = find(stdl<minstd); %only look for acceptable std

[~,J] = min(abs(zs(idx))); % take the min zscore

y = seq(idx(J),:);

evol = sum(seq,2)/ntrl;

    function s = seqswitch(v,csrt,ntrl)
        s = zeros(1,ntrl);
        for t = csrt+1:ntrl
            s(t) =  (rand(1)<=v);
            if (s(t) == 1) && sum(s((t-csrt):(t-1))~=0)
                s(t) = 0;
            end
        end
    end

end