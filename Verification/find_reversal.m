function block_rev = find_reversal(expe,vol)
% function giving  reversal index for a given blck, volatility level in expe
% assuming same number of trials per volatility level
% removing inter-vol reversals

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



