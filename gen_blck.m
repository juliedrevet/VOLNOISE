function [blck] = gen_blck(volperm)
%  GEN_BLCK  Generate block of VOLNOISE experiment
%
%  Usage: [blck] = GEN_BLCK(volperm)
%
%  Output: BLCK structure containing following fields
%   * volperm      - volatility levels permutation => 1:low | 2:middle | 3: high volatility
%
%   NEW FIELDS:
%   * ntrl       - number of trial for this block
%   * switch_seq - reward prob switch => 1: switch or 0: no switch
%   * switch_spe - same as switch_seq but with 2: intervolatility switch
%   * false_seq  - false positives sequence => 1: false positive or 0: no
%   * reward_seq - probability of reward at each trial => array(1,ntrl)
%   * outcome    - actual rewarding sequence => 1: shape(1) or 2: shape(2)

if nargin < 1
    error('Missing configuration structure!');
end

% create block output structure
blck         = [];
blck.volperm = volperm;
blck.ntrl    = 0;

last_epi  = zeros(1,length(volperm)-1);
first_epi = zeros(1,length(volperm)-1);

for ivol = 1:length(volperm)
    if volperm(ivol) == 0
        fprintf('  * generating episodes for training volatility...\n');
    else
        fprintf('  * generating episodes for volatility level %d...\n',volperm(ivol));
    end
    % set subblock parameters
    if volperm(ivol) == 0 % practice block >> OK?
        davg     = 20;
        nepi     = 3;
        p_reward = .9;
    elseif volperm(ivol) == 1 % low volatility
        davg     = 20;
        nepi     = 3;
        p_reward = .85;
    elseif volperm(ivol) == 2 % middle volatility
        davg     = 15;
        nepi     = 4;
        p_reward = .85;        
    elseif volperm(ivol) == 3 % high volatility
        davg     = 10;
        nepi     = 6;
        p_reward = .85; 
    end
           
    ntrl   = davg*nepi;  % number of trials per volatility level
    dlim   = [5 ntrl];   % min|max number of trials before reversal
    rfalse = 1-p_reward; % proportion of false positives

    % generate subblock episodes (Valentin's function)
    b = gen_epi(davg,dlim,nepi);   
    
    blck.nepi(ivol) = nepi;
    blck.davg(ivol) = davg;
    blck.dlim       = dlim;
    blck.prev(ivol) = b.prev;
    blck.pr(ivol)   = b.pr;
    blck.pb(ivol)   = b.pb;    
    blck.correct(ivol,:) = b.ys;

    switch_idx = cumsum(b.xs)+1; % idx of switch
    switch_seq = zeros(1,ntrl);
    switch_seq(switch_idx(1:end-1)) = 1;
    blck.switch_seq(ivol,:) = switch_seq;
    blck.switch_spe(ivol,:) = switch_seq;
    
    if ivol < length(volperm)
        last_epi(ivol) = b.xs(end);
    end
    if ivol > 1
        first_epi(ivol-1) = b.xs(1);
    end
  
    blck.ntrl = blck.ntrl+ntrl;
end

% combining subblocks
fprintf('  * combining block structure...\n');
% generate inter-volatility reversal (last 2 : first 2) and check for sequence continuity
if length(volperm) > 1
    for i = 1:(length(volperm)-1)
        interrnd = 2*rand-1;
        interrev =round(interrnd);
        while (last_epi(i)+interrev < 5)||(first_epi(i)-interrev < 5)
            interrnd = 2*rand-1;
            interrev =round(interrnd);
        end
        if interrnd < 0
            blck.switch_seq(i,end+interrev) = 1;
            blck.switch_spe(i,end+interrev) = 2;
            blck.correct(i,(end+interrev):end)   = 3-blck.correct(i,(end+interrev):end);
            if blck.correct(i+1,1) ~= blck.correct(i,end)
                % reverse all following sequences
                for k = (i+1):length(volperm)
                    blck.correct(k,:) = 3-blck.correct(k,:);
                end
            end
        elseif interrnd > 0
            if blck.correct(i,end) == blck.correct(i+1,1)
                % reverse all following sequences
                for k = (i+1):length(volperm)
                    blck.correct(k,:) = 3-blck.correct(k,:);
                end
            end
            blck.switch_seq(i+1,1+interrev) = 1;
            blck.switch_spe(i+1,1+interrev) = 2;
            blck.correct(i+1,1:interrev)    = 3-blck.correct(i+1,1:interrev);
        end
    end
end

% generating false positive sequences
for i = 1:length(volperm)
    switch_seq = blck.switch_seq(i,:);
    %false positive sequence {0;1}
    false_seq = sample_false; % same number of false positives at 1st and 2nd half
    idx_false = find(false_seq);
    % NOT: more than 2 consecutives false positives, false positive at the beginning
    % or at the end of the block, just before or on top of a switch
    while any(diff(idx_false(diff(idx_false)==1))==1)||ismember(1,idx_false)||ismember(ntrl,idx_false)... % maybe remove following conditions
            ||any(switch_seq(idx_false))||any(switch_seq(idx_false+1))
        false_seq = sample_false;
        idx_false = find(false_seq);
    end
    blck.false_seq(i,:) = false_seq;
end

blck.switch_seq = reshape(blck.switch_seq',1,length(volperm)*ntrl);
blck.switch_spe = reshape(blck.switch_spe',1,length(volperm)*ntrl);
blck.false_seq  = reshape(blck.false_seq',1,length(volperm)*ntrl);
blck.correct    = reshape(blck.correct',1,length(volperm)*ntrl);
blck.reward_seq = (blck.correct==1)*p_reward+(blck.correct==2)*(1-p_reward);
blck.outcome    = blck.correct;
blck.outcome(logical(blck.false_seq)) = 3-blck.outcome(logical(blck.false_seq));

    function s = sample_false
        if ~mod(int8(rfalse*ntrl),2)
            s = [Shuffle(kron([1 zeros(1,int8(1/rfalse-1))],ones(1,int8(rfalse*ntrl/2))))...
                Shuffle(kron([1 zeros(1,int8(1/rfalse-1))],ones(1,int8(rfalse*ntrl/2))))];
        else % ntrl*rfalse odd
            s = Shuffle([ones(1,int8(ntrl*rfalse)) zeros(1,ntrl*(1-rfalse))]);
            while ((sum(s(1:(ntrl/2)))+1)~=sum(s((ntrl/2)+1:end)))&&(sum(s(1:(ntrl/2)))~=(sum(s((ntrl/2)+1:end))+1))
                s = Shuffle([ones(1,int8(ntrl*rfalse)) zeros(1,ntrl*(1-rfalse))]);
            end
        end
    end

end