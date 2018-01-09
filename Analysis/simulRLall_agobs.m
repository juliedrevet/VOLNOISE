function simulRLall_agobs(subj, n, fun)

if nargin<3
    n = 1000;
    fun = @fminsearch;
elseif nargin == 3
    fun = @fminsearch;
end

subj_data = load_data(~ismember(1:23,subj));

for blck = 3:10
    expe = subj_data.expe(blck);
    
    taskid = expe.cfg.taskid(blck);
    condtn = expe.cfg.condtn(blck);
    
    [alphabeta_olow, alphabeta_alow, alphabeta_ohigh, alphabeta_ahigh] = findab4(subj,fun);
    
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
    
    fprintf('Block %d alpha: %5.2f\n', blck-2, alpha)
    fprintf('Block %d beta: %5.2f\n', blck-2, beta)
    
    out = deltaQlearn(alpha,outcome);
    
    actionsbl = (rand(length(out),n)<(1./(1+exp(beta*out))))*1.;
    actionsbl(actionsbl == 0) = 2;
    
    actions((blck-3)*120+1:(blck-2)*120,:) = actionsbl;
    resp((blck-3)*120+1:(blck-2)*120) = expe.rslt.resp;
    correct((blck-3)*120+1:(blck-2)*120) = expe.rslt.correct;
    false_seq((blck-3)*120+1:(blck-2)*120) = expe.blck.false_seq;
end

figure();
plot(mean(actions,2))
hold on
plot(resp,'LineWidth',1.5)
plot(correct,'--','LineWidth',2)
ylim([.8 2.2]);
yticks([1 2]);
plot(false_seq*.9,'*')
legend({'RL model','subject','correct','false positives'})
xlim([1 960]);

end