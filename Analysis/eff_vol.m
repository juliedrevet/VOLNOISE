function eff_vol(v,csrt,nsubj)


% v is an array of volatilities: 0.0100    0.0172    0.0295    0.0508    0.0873    0.1500
n = 5000;
ntrl = 120;
chosen_seq = zeros(length(v),ntrl,nsubj);
all_evol = zeros(length(v),n,nsubj);

k = 1;
for vol = v
    for s = 1:nsubj
        [y, evol] = gen_switch(n,ntrl,vol,csrt,30);
        chosen_seq(k,:,s) = y;
        all_evol(k,:,s)   = evol;
        
    end
    k = k+1;
end

all_evol = reshape(all_evol,length(v),nsubj*n);

%% plot
% plot mean length before switch and std boundaries for each volatility
z = mean(all_evol,2); % y axis
zu = z+std(all_evol,0,2); % upper bound
zl = z-std(all_evol,0,2); % lower bound
xvol = v; % x axis

figure()
hold all

%shaded std area of our generated samples
p = fill([xvol flip(xvol)], [zu' flip(zl)'], 'blue', 'linestyle', 'none');
set(p,'facealpha',.1) % change opacity
set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % don't show in legend
plot(xvol,z,'blue','LineWidth',.5)
plot(xvol,mean(sum(chosen_seq,2),3)/ntrl,'blue','LineWidth',2)
plot(xvol,mean(sum(chosen_seq,2),3)/ntrl,'b+')
plot(xvol,xvol,'ro')

xlim([min(xvol) max(xvol)]);
xlabel('volatility','FontSize',14);
ylabel('effective volatility  ( n_{reversals}/ n_{trials} )','FontSize',14);
title(sprintf('Effective volatility with minimum %d trials before switch',csrt))

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
