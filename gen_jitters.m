function jitter_seq = gen_jitters(subj,session,expe)
% generate jitter sequence for subj

% jitter values in seconds
%iti   = [.5 2.2 3.9]; % inter-trial interval
%stout = [.4 2.1 3.8]; % stimulus end - outcome interval

% 9 possible combinations
%          ITI STOUT
jitters = [1    1; ... % 1
           1    2; ... % 2
           1    3; ... % 3
           2    1; ... % 4
           2    2; ... % 5
           2    3; ... % 6
           3    1; ... % 7
           3    2; ... % 8
           3    3];    % 9

lsq = latin_square; % dimension 3 (3 blocks, 3 trials before|after reversal)

%% take parameters from expe
nblck = length(expe.blck);
ntrl  = expe.blck(1).ntrl;
nvol  = length(expe.blck(1).volperm);
ntrl_vol = ntrl/nvol;

%% generate jitters for beginning of block >> control on first 6 trials
iti   = zeros(3,6);
stout = zeros(3,6);

jitter_seq.iti   = zeros(nblck,ntrl);
jitter_seq.stout = zeros(size(jitter_seq.iti));

if session == 1
    tempsubj = subj;
else %session 2
    tempsubj = subj+4;
end

iti(:,1:3)   = lsq(:,:,mod(tempsubj-1,12)+1);
stout(:,1:3) = lsq(:,:,mod(tempsubj+0,12)+1);
iti(:,4:6)   = lsq(:,:,mod(tempsubj+1,12)+1);
stout(:,4:6) = lsq(:,:,mod(tempsubj+2,12)+1);

% fix beginning jitters
jitter_seq.iti(:,1:6)   = iti;
jitter_seq.stout(:,1:6) = stout;

%% generate jitters before|after reversal per volatility levels (all 3 blocks at a time through latin square)
for ivol = 1:nvol % for each volatility level
    idx_vol = find_reversal(expe,ivol);
    lsq_vol = size(idx_vol,2); % number of latin squares needed
    combisubjpre  = genjitters(tempsubj,lsq_vol); % generate pre-reversal 1:9 combi
    combisubjpost = genjitters(tempsubj,lsq_vol); % generate post-reversal 1:9 combi
    for l = 1:lsq_vol % for each reversal
        iti_pre   = itistout2iti(combisubjpre(:,:,l));
        iti_post  = itistout2iti(combisubjpost(:,:,l));
        stout_pre  = itistout2stout(combisubjpre(:,:,l));
        stout_post = itistout2stout(combisubjpost(:,:,l));
        for iblck = 1:nblck % assign iti and stout values per block
            jitter_seq.iti(iblck,(idx_vol(iblck,l)-3):(idx_vol(iblck,l)-1)) = iti_pre(iblck,:);
            jitter_seq.iti(iblck,(idx_vol(iblck,l)):(idx_vol(iblck,l)+2))   = iti_post(iblck,:);
            jitter_seq.stout(iblck,(idx_vol(iblck,l)-3):(idx_vol(iblck,l)-1)) = stout_pre(iblck,:);
            jitter_seq.stout(iblck,(idx_vol(iblck,l)):(idx_vol(iblck,l)+2))   = stout_post(iblck,:);
        end
    end
end

%% generate jitters for remaining trials
% control: have uniformed distribution of 1:9 combinations per volatility level

for iblck = 1:nblck
    for k = 1:nvol % per volatility
        idx_empty = find(jitter_seq.iti(iblck,((k-1)*ntrl_vol+1):(k*ntrl_vol)) == 0);
        idx_full  = find(jitter_seq.iti(iblck,((k-1)*ntrl_vol+1):(k*ntrl_vol)) ~= 0);
        is = give_combination(jitter_seq.iti(iblck,idx_full+(k-1)*ntrl_vol),jitter_seq.stout(iblck,idx_full+(k-1)*ntrl_vol));
        h = hist(is,1:9);
        jitts = kron(1:9,ones(ceil(ntrl_vol/9),1));
        for l = 1:9
            jitts(1:h(l),l) = 0;
        end
        jitts = nonzeros(jitts)';
        jitts = Shuffle(jitts);
        % not 3 consecutives
        while any(diff(abs(diff(jitts)))==0)
            jitts = Shuffle(jitts);
        end
        jitter_seq.iti(iblck,idx_empty+(k-1)*ntrl_vol)   = itistout2iti(jitts(1:length(idx_empty)));
        jitter_seq.stout(iblck,idx_empty+(k-1)*ntrl_vol) = itistout2stout(jitts(1:length(idx_empty)));
    end   
end
    
    function iti = itistout2iti(itistout)
        iti = zeros(size(itistout));
        for i = 1:size(itistout,1)
            for j = 1:size(itistout,2)
                iti(i,j) = jitters(itistout(i,j),1);
            end
        end
    end

    function stout = itistout2stout(itistout)
        stout = zeros(size(itistout));
        for i = 1:size(itistout,1)
            for j = 1:size(itistout,2)
                stout(i,j) = jitters(itistout(i,j),2);
            end
        end
    end

end

function block_rev = find_reversal(expe,vol)
% function giving  reversal index per blck for a given volatility level in expe
% assuming same number of trials per volatility level
% removing inter-volatility-levels reversals

for iblck = 1:length(expe.blck)
    blck     = expe.blck(iblck);
    ntrl     = blck.ntrl;
    nvol     = length(blck.volperm);
    ntrl_vol = ntrl/nvol;
    ivol     = find(blck.volperm == vol);
    irev     = find(blck.switch_seq(((ivol-1)*ntrl_vol+1):ivol*ntrl_vol));
    block_rev(iblck,:) = irev(ismember(irev,3:ntrl_vol-2))+(ivol-1)*ntrl_vol; % remove inter-reversals
end

end

function itistout = give_combination(iti,stout)
% return iti-stout jitters combinations (1:9) given iti and stout arrays
% (size not necessarily 3*3)

% 9 possible combinations
%          ITI STOUT
combi = [11; ... %1
         12; ... %2
         13; ... %3
         21; ... %4
         22; ... %5
         23; ... %6
         31; ... %7
         32; ... %8
         33];    %9
c = cell(size(iti));
itistout = zeros(size(iti));

for j = 1:size(iti,1)
    for k = 1:size(iti,2)
        c{j,k}=sprintf('%d%d',iti(j,k),stout(j,k));
        itistout(j,k) = find(str2double(c{j,k})==combi);
    end
end

end