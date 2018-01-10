function expe = orga_volnoise(ex)
% modification essai pour github


%ntrl = length(ex.blck)*ex.blck(1).ntrl;
expe = [];

itrl_session = 1;
for iblck = 1:length(ex.blck)
    for itrl = 1:ex.blck(iblck).ntrl
        expe(itrl_session).subj    = ex.hdr.subj;
        expe(itrl_session).session = ex.hdr.session;
        expe(itrl_session).trial   = itrl_session;
        expe(itrl_session).blck    = iblck;
        expe(itrl_session).trial_blck = itrl;
        expe(itrl_session).volatility = ex.blck(iblck).volperm(ceil(itrl/60));
        expe(itrl_session).stim1 = ex.stim(iblck).shape(1);
        expe(itrl_session).stim2 = ex.stim(iblck).shape(2);
        expe(itrl_session).exp_resp = ex.blck(iblck).correct(itrl);
        expe(itrl_session).exp_shape = ex.stim(iblck).shape(ex.blck(iblck).correct(itrl));
        expe(itrl_session).trap_trial = ex.blck(iblck).false_seq(itrl);
        expe(itrl_session).switch_trial = ex.blck(iblck).switch_seq(itrl);
        expe(itrl_session).jitter_stimfb = ex.clck(iblck).toutc(itrl)-ex.clck(iblck).tstim(itrl);
        if itrl ~= ex.blck(iblck).ntrl
            expe(itrl_session).jitter_iti = ex.clck(iblck).tstim(itrl+1)-ex.clck(iblck).toutc(itrl);
        end
        expe(itrl_session).resp = ex.rslt(iblck).resp(itrl);
        expe(itrl_session).resp_key = ex.rslt(iblck).respkb(itrl);
        expe(itrl_session).correct = double(ex.blck(iblck).outcome(itrl) == ex.rslt(iblck).resp(itrl));
        if ex.rslt(iblck).resp(itrl)==1
            expe(itrl_session).reward = ex.blck(iblck).reward_seq(itrl);
        elseif ex.rslt(iblck).resp(itrl)==2
            expe(itrl_session).reward = 1-ex.blck(iblck).reward_seq(itrl);
        end
        expe(itrl_session).rt = ex.rslt(iblck).rt(itrl);
        expe(itrl_session).fb = ex.stim(iblck).feedback(itrl);
        
        itrl_session = itrl_session+1;
    end

end

end