function verif_jitters(expe)

nblck = length(expe.stim);
ntrl = length(expe.stim(1));

iti   = [.5 2.2 3.9]; % inter-trial interval (fixation point before stim onset)
stout = [.4 2.1 3.8]-.3; % end of stim to outcome duration

% all jitters equally distributed per volatility level
for iblck = 1:nblck
    for ivol = 1:length(expe.blck(iblck).volperm)
        idx_vol = expe.blck(iblck).volperm(ivol);
        portion = ((ivol-1)*60+1) : ivol*60;
        stouts(iblck,:,idx_vol) = expe.clck(iblck).toutc(portion)-expe.clck(iblck).tstim(portion)-1.3;
        itis(iblck,:,idx_vol) = expe.clck(iblck).tstim(portion)-expe.clck(iblck).tfix1(portion);
% % %         itistouts(iblck,:,idx_vol) = give_combi(expe.stim(iblck).iti_jitters_idx(portion),...
% % %             expe.stim(iblck).stout_jitters_idx(portion));
    end
end

h_iti = hist(itis,iti);
h_stout = hist(stout,stout);
% % % h_combi = hist(itistouts,1:9);







end


% % % for k = 1:1%25
% % %     for session = 1:2
% % %         expe = gen_expe(k,session);
% % %         
% % %         datapath = sprintf('./Data/');
% % %         filename = sprintf('VOLNOISE_IRM_expe%02d_session%02d',k,session);
% % %         filename = [filename,'.mat'];
% % %         filename = fullfile(datapath,filename);
% % %         save(filename,'expe');
% % %         
% % %         nvol  = length(expe.blck(1).volperm);
% % %         nblck = length(expe.blck);
% % %         jitter_seq = gen_jitters(k,session,expe);
% % %         
% % %         iti = [];
% % %         stout = [];
% % %         
% % %         for ivol = 1:nvol
% % %             idx_vol = find_reversal(expe,ivol);
% % %             iti_pre  = [];
% % %             iti_post = [];
% % %             stout_pre  = [];
% % %             stout_post = [];
% % %             for iblck = 1:nblck
% % %                 for j = 1:size(idx_vol,2)
% % %                     iti_pre = [iti_pre; jitter_seq.iti(iblck,(idx_vol(iblck,j)-3):(idx_vol(iblck,j)-1))];
% % %                     iti_post = [iti_post; jitter_seq.iti(iblck,(idx_vol(iblck,j)):(idx_vol(iblck,j)+2))];
% % %                     stout_pre = [stout_pre; jitter_seq.stout(iblck,(idx_vol(iblck,j)-3):(idx_vol(iblck,j)-1))];
% % %                     stout_post = [stout_post; jitter_seq.stout(iblck,(idx_vol(iblck,j)):(idx_vol(iblck,j)+2))];
% % %                 end
% % %             end
% % %             iti(ivol).pre = iti_pre;
% % %             iti(ivol).post = iti_post;
% % %             stout(ivol).pre = stout_pre;
% % %             stout(ivol).post = stout_post;
% % %         end
% % %         filename = sprintf('VOLNOISE_IRM_itiS%02d_session%02d',k,session);
% % %         filename = [filename,'.mat'];
% % %         filename = fullfile(datapath,filename);
% % %         save(filename,'iti');
% % %         filename = sprintf('VOLNOISE_IRM_stoutS%02d_session%02d',k,session);
% % %         filename = [filename,'.mat'];
% % %         filename = fullfile(datapath,filename);
% % %         save(filename,'stout');
% % %     end
% % % end