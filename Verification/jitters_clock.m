function [passed_iti, passed_stout, match_iti,match_stout] = jitters_clock(expe,epsilon)
% this function verifies that the actual flipping time corresponds to the
% desired jitter duration for ITI (inter-trial-interval) and STOUT (stim to 
% outcome) jitters with a time tolerance of epsilon seconds (default 1e-2)
% it returns 2 logical arrays containing 1 if flip done with the correct timing
% match_iti and match_stout (nblck * ntrl)


if nargin<1
    error('Missing experiment structure!');
elseif nargin == 1
    epsilon = expe.stim(1).ifi/2;
end    

ntrl  = length(expe.rslt(1).resp);
nblck = length(expe.blck);
match_iti   = zeros(nblck,ntrl);
match_stout = zeros(nblck,ntrl);

for iblck = 1:nblck
    for itrl = 1:ntrl
        lr_iti = expe.stim(iblck).iti_jitters(itrl)-epsilon; % lower range
        ur_iti = expe.stim(iblck).iti_jitters(itrl)+epsilon; % upper range
        lr_stout = expe.stim(iblck).stout_jitters(itrl)-epsilon;
        ur_stout = expe.stim(iblck).stout_jitters(itrl)+epsilon;
        iti   = expe.clck(iblck).tstim(itrl)-expe.clck(iblck).tfix1(itrl);
        stout = expe.clck(iblck).toutc(itrl)-expe.clck(iblck).tstim(itrl)-1.3;
        match_iti(iblck,itrl) = (iti>=lr_iti)&(iti<=ur_iti);
        match_stout(iblck,itrl) = (stout>=lr_stout)&(stout<=ur_stout);
    end
end
passed_iti   = sum(match_iti(:))==nblck*ntrl;
passed_stout = sum(match_stout(:))==nblck*ntrl;

end
        