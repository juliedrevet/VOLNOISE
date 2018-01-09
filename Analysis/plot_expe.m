function plot_expe(expe)

colors = [0      0.4470 0.7410;... % blue (dark)
    0.8500 0.3250 0.0980;... % red (light)
    0.9290 0.6940 0.1250;... % yellow
    0.4940 0.1840 0.5560;... % purple
    0.4660 0.6740 0.1880;... % green
    0.3010 0.7450 0.9330;... % blue (light)
    0.6350 0.0780 0.1840;... % red (dark)
    ];



nblck = length(expe.blck);
ntrl  = expe.blck(1).ntrl;
nvol  = length(expe.blck(1).volperm);
ntrl_vol = ntrl/nvol;

colors_vol = parula;
%colors_vol = colors_vol(linspace(5,10*nvol,nvol),:);
colors_vol = colors_vol(ceil(linspace(10,25,nvol)),:);

h = figure();

for iblck = 1:nblck
    hh(iblck) = subplot(nblck,1,iblck);
    blck = expe.blck(iblck);
    p_reward = [max(blck.reward_seq) min(blck.reward_seq)];
    p1=plot(1:ntrl,blck.reward_seq,'Color',colors_vol(1,:),'Linewidth',2);
    set(get(get(p1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % don't show in legend
    for ivol = double(blck.volperm)
        hold on
        xvol = ((ivol-1)*ntrl_vol+1):ivol*ntrl_vol;
        yrewmin = (min(p_reward)-.1)*ones(1,ntrl_vol);
        yrewmax = (max(p_reward)+.1)*ones(1,ntrl_vol);
        p = fill([xvol fliplr(xvol)], [yrewmin fliplr(yrewmax)],colors_vol(blck.volperm(ivol),:),'linestyle','none');
        %text(min(xvol)+25,min(volatilities)*0.5, 'Obs', 'FontSize',10)
        set(p,'facealpha',.08)
        set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
    
    plot(find(blck.false_seq),.95*ones(sum(blck.false_seq)),'Color',colors(2,:),'Marker','.','LineStyle','none')
    ylim([.1,1]);   yticks(sort(p_reward));
    yticklabels({sprintf('%6.2f',p_reward(2)),sprintf('%6.2f',p_reward(1))});
    xlim([0,ntrl]); xticks(0:ntrl_vol:ntrl); xticklabels('');
    
    if iblck == 2
        ylabel('reward probability','FontSize',14)
    elseif iblck == 3
% % %         l1 = legend('false positive','Location','southeast');
% % %         l1.Box='off';
% % %         l1.FontSize = 14;
        xt = 0:ntrl_vol:ntrl;
        xlim([0,ntrl]); xticks(xt); xlabel('trial','FontSize',14);
        xtl = cell(1,length(xt));
        for i = 1:length(xt)
            xtl{i} = sprintf('%d',xt(i));
        end
        xticklabels(xtl);
        
    end

%h1.Units = 'pixels';
%h2.Units = 'pixels';
h_in(iblck,:) = hh(iblck).Position; % [left bottom width height]
h_out(iblck,:) = hh(iblck).OuterPosition;

    
end

% bottom pos
bottpos = h_in(3,2);
% top pos
toppos  = h_in(1,2)+h_in(1,4);

newheight = (toppos-bottpos)/3;

hh(1).Position = [h_in(1,1) toppos-newheight h_in(1,3) newheight];
hh(2).Position = [h_in(2,1) toppos-2*newheight h_in(2,3) newheight];
hh(3).Position = [h_in(3,1) h_in(3,2) h_in(3,3) newheight];
end