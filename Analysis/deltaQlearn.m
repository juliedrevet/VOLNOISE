function out = deltaQlearn(alpha,outcome,rst,reward)

N=length(outcome); %number of trials

if nargin<3
    rst = N;
    reward = 1;
elseif nargin == 3
    reward = 1;
end

deltaQ = zeros(N,1); % = Q2-Q1
deltaQ(1) = 0;

for ind = 1:(N - 1)
    if mod(ind,rst)==0
        deltaQ(ind+1) = 0; %reset deltaQ
    else
        if ((outcome(ind) == 2)) % delta_R = +reward
            deltaQ(ind+1) = deltaQ(ind) + alpha*(reward-deltaQ(ind));
        else  % delta_R = -reward
            deltaQ(ind+1) = deltaQ(ind) + alpha*(-reward-deltaQ(ind));
        end
    end
end

out = deltaQ;

end