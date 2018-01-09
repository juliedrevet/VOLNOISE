function jitter = verif(expe,jitter_seq,k)

if nargin == 2
    k = 3;
elseif (nargin == 3)
    if (abs(k)>5)
        error(' Only possible to look %d trials back',expe.blck(1).dlim(1));
    elseif k == 0
        error('0 is not an acceptable argument');
    end
end

nvol  = length(expe.blck(1).volperm);
nblck = length(expe.blck);

%jitter = struct(nvol);

for ivol = 1:nvol
    idx_vol = find_reversal(expe,ivol);
kiti  = [];
kstout = [];
    for iblck = 1:nblck
        for j = 1:size(idx_vol,2)
            if k<0
                kiti = [kiti;...
                    jitter_seq.iti(iblck,(idx_vol(iblck,j)+k):(idx_vol(iblck,j)-1))];
                kstout = [kstout;...
                    jitter_seq.stout(iblck,(idx_vol(iblck,j)+k):(idx_vol(iblck,j)-1))];

            else  
                kiti = [kiti; jitter_seq.iti(iblck,(idx_vol(iblck,j)):(idx_vol(iblck,j)+k-1))];
                kstout = [kstout; jitter_seq.stout(iblck,(idx_vol(iblck,j)):(idx_vol(iblck,j)+k-1))];
            end
        end
    end
    
    jitter(ivol).iti   = kiti;
    jitter(ivol).stout = kstout;
    h(:,:,ivol)=hist(give_combi(jitter(ivol).iti,jitter(ivol).stout),1:9);
end

end