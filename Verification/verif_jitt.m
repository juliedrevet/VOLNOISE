function verif_jitt(expe)

iti_val   = [.5 2.2 3.9]; % inter-trial interval (fixation point before stim onset)
stout_val = [.4 2.1 3.8]-.3;


for iblck = 1:length(expe.blck)
    iti(iblck,:)   = expe.stim(iblck).iti_jitters;
    stout(iblck,:) = expe.stim(iblck).stout_jitters;
    iti(iblck,iti(iblck,:)<iti_val(2)) = 1;
    iti(iblck,iti(iblck,:)>iti_val(2)) = 3;
    iti(iblck,iti(iblck,:)==iti_val(2)) = 2;
    stout(iblck,stout(iblck,:)<stout_val(2)) = 1;
    stout(iblck,stout(iblck,:)>stout_val(2)) = 3;
    stout(iblck,stout(iblck,:)==stout_val(2)) = 2;
end

jitter_seq.iti = iti;
jitter_seq.stout = stout;

verif(expe,jitter_seq,-3)
verif(expe,jitter_seq,3)


end
