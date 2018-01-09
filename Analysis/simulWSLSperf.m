function perfWSLS = simulWSLSperf(blck,n)

if nargin==1
    n = blck.ntrl; 
end


% win-stay lose-switch:
wsls_ans    = zeros(1,n); %{1;2}
wsls_ans(1) = randi(2);
for i = 2:n
    if wsls_ans(i-1) == blck.outcome(i-1)
        wsls_ans(i) = wsls_ans(i-1); 
    else
        wsls_ans(i) = 3-wsls_ans(i-1);
    end
end

perfWSLS = sum(wsls_ans==blck.outcome(1:n))/n;

end