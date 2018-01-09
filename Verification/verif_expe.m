function verif_expe(expe)

addpath ../Analysis
addpath ../Toolboxes/Rand

plot_expe(expe)

% check volatility levels permutation (latin square?)
vol_perm = zeros(3,3);
for iblck = 1:3
    vol_perm(iblck,:) = expe.blck(iblck).volperm;
end
vol_perm
fprintf('\n');

% check amount of false positives
for iblck = 1:3
    for ivol = 1:3
        fprintf('Block %d, volatility level %d: %d false positives\n',...
            iblck, expe.blck(iblck).volperm(ivol), sum(expe.blck(iblck).false_seq((ivol-1)*60+1:ivol*60)));
    end
end
fprintf('\n')

% check not more than 2 consecutive false positives
for iblck = 1:3
    [xcon,ncon] = GetConsecutiveValues(expe.blck(iblck).false_seq);
    fprintf('Block %d: maximal %d consecutive false positives\n',iblck,max(ncon(xcon==1)));
end
fprintf('\n')

% not false feedback at the trial before or at a reversal
for iblck = 1:3
    idx_rev = find(expe.blck(iblck).switch_seq);
    fprintf('Block %d: number of false positives on the trial before or at a reversal: %d\n',...
        iblck,sum(expe.blck(iblck).false_seq(idx_rev))+sum(expe.blck(iblck).false_seq(idx_rev-1)))
end
fprintf('\n')

% check number of reversals per bloc: max 12
for iblck = 1:3
    if sum(expe.blck(iblck).switch_seq)== 12
        fprintf('Block %d: %d reversals ... OK\n',iblck, sum(expe.blck(iblck).switch_seq));
    else
        fprintf('Block %d: %d reversals ... not correct\n',iblck, sum(expe.blck(iblck).switch_seq));
    end
end
fprintf('\n');

% % % % check number of episodes per volatility level
% % % for iblck = 1:3
% % %     for ivol = 1:3
% % %         if ivol == 1
% % %             portion = 1:62;
% % %         elseif ivol == 2
% % %             portion = 59:122;
% % %         else
% % %             portion = 119:180;
% % %         end   
% % %         fprintf('Block %d, volatility level %d: %d episodes\n',...
% % %             iblck, expe.blck(iblck).volperm(ivol), numel(find([expe.blck(iblck).switch_seq(portion) 0])>=5));
% % %         diff(find(expe.blck(iblck).switch_seq(portion)))
% % %         sum(find([expe.blck(iblck).switch_seq(portion) 0])>=5)
% % %         %find(expe.blck(iblck).switch_seq(portion))
% % %     end
% % % end
% % % fprintf('\n')  

% check length of episodes: not smaller than 5 trials
for iblck = 1:3
    fprintf('Block %d: minimum episode length is %d\n',iblck,min(diff([1 find(expe.blck(iblck).switch_seq) 181])));
end
fprintf('\n')

% check stimulus: shape not more than 3 times on the same side
for iblck = 1:3
    [~,ncon] = GetConsecutiveValues(expe.stim(iblck).pos(1,:));
    fprintf('Block %d: A given shape is maximally %d consecutive trials on the same left (respectively right) side\n',iblck, max(ncon));
end