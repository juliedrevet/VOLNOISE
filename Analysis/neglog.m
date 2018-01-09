function nllh = neglog(Q,data,beta)

N = size(Q,1);
if length(data) ~= N
    error('vectors must be the same length please check!')
else
    
    idx1 = find(data == 1);
    idx2 = find(data == 2);
    
    if size(Q,2) ~=1
        deltaQ = Q(:,2)-Q(:,1);
    else
        deltaQ = Q;
    end
    l(idx1) = 1./(1+exp(beta*deltaQ(idx1)));
    l(idx2) = 1./(1+exp(-beta*deltaQ(idx2)));
    nllh = (-1)*sum(log(l(l~=0))); % remove when no choice made
end

end
