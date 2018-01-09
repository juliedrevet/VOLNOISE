function simulRL_agobs(subj, blck, n, fun)

if nargin<3
    n = 1000;
    fun = @fmincon;
elseif nargin == 3
    fun = @fmincon;
end

subj_data = load_data(~ismember(1:23,subj));
expe = subj_data.expe(blck);

taskid = expe.cfg.taskid(blck);
condtn = expe.cfg.condtn(blck);

[alphabeta_olow, alphabeta_alow, alphabeta_ohigh, alphabeta_ahigh] = findab4(subj,fun);

[alphabeta,~,~] = findab(subj,fun);
alpha2 = alphabeta(blck,1);
beta2 = alphabeta(blck,2);

if taskid==1 % not dependent on choice
    if condtn == 1
        alpha = alphabeta_olow(1);
        beta  = alphabeta_olow(2);
    else
        alpha = alphabeta_ohigh(1);
        beta  = alphabeta_ohigh(2);
    end
    outcome = expe.blck.outcome(1,:);
elseif taskid==2 % dependent on choice
    if condtn == 1
        alpha = alphabeta_alow(1);
        beta  = alphabeta_alow(2);
    else
        alpha = alphabeta_ahigh(1);
        beta  = alphabeta_ahigh(2);
    end
    outcome = expe.blck.epimap*expe.blck.color_seq...
        +(3-expe.blck.epimap)*(~expe.blck.color_seq);
end

fprintf('alpha1: %5.2f\n', alpha)
fprintf('beta1: %5.2f\n', beta)
fprintf('alpha2: %5.2f\n', alpha2)
fprintf('beta2: %5.2f\n', beta2)

out = deltaQlearn(alpha,outcome);
out2 = deltaQlearn(alpha2,outcome);

actions = (rand(length(out),n)<(1./(1+exp(beta*out))))*1.;
actions(actions == 0) = 2;

actions2 = (rand(length(out2),n)<(1./(1+exp(beta2*out2))))*1.;
actions2(actions2 == 0) = 2;

% % % %mean neg likelihod of simulated actions
% % % neglog(deltaQlearn(alpha,expe.rslt.resp,120),expe.rslt.resp,beta)
% % % 
% % % ll = zeros(1,n);
% % % for i = 1:n
% % %     ll(i) = neglog(deltaQlearn(alpha,expe.rslt.resp,120),actions(:,i)',beta);
% % % end
% % % 
% % % mean(ll)

perf = expe.rslt.perf;
sim_perf = sum((actions' == expe.rslt.correct),2)/120;
sim_perf2 = sum((actions2' == expe.rslt.correct),2)/120;

figure('Color','white');
set(gca,'Layer','top','Box','off');
plot(mean(actions,2))
hold on
plot(expe.rslt.resp,'LineWidth',2)
plot(expe.rslt.correct,'--','LineWidth',2)
ylim([.8 2.2]);
yticks([1 2]);
plot(expe.blck.false_seq*.9,'*')
set(gca,'Layer','top','Box','off');
set(gca,'Visible','off')
legend({'RL model','subject','correct','false positives'})
legend('boxoff')
set(gca,'FontName','Helvetica','FontSize',16); 


figure()
subplot(2,1,1)
histogram(sim_perf,'Normalization','probability')
xlim([min(min(sim_perf),min(sim_perf2)) max(max(sim_perf),max(sim_perf2))]);
hold on
plot(perf,0,'ro')

subplot(2,1,2)
histogram(sim_perf2,'Normalization','probability')
hold on
plot(perf,0,'ro')
xlim([min(min(sim_perf),min(sim_perf2)) max(max(sim_perf),max(sim_perf2))]);

end