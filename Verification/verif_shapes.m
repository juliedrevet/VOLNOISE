function [passed_all, passed_side, passed_eq, passed_eq_rew] = verif_shapes(expe)
% verification regarding shapes presentation

nblck = length(expe.stim);
ntrl = length(expe.stim(1));

% all shapes presented once (all shapes different)
allshapes = zeros(nblck,2);
for iblck=1:nblck
     % just look at left presented shape it's enough
    allshapes(iblck,:) = unique(expe.stim(iblck).pos_shape(1,:));
end
passed_all = (length(unique(allshapes)) == nblck*2);

% a given shape only presented max 3 times on the same side
consec_side = zeros(1,nblck);
for iblck = 1:nblck
    [~,ncon] = GetConsecutiveValues(expe.stim(iblck).pos_shape(1,:));
    consec_side(iblck) = max(ncon);
end
passed_side = all((consec_side <= 3*ones(1,nblck)));

% both shapes equally presented on the left and the right side
h = zeros(nblck,2);
for iblck = 1:nblck
    h(iblck,:) = hist(expe.stim(iblck).pos_shape(1,:),2);
end
passed_eq = all(h(:,1)==h(:,2));

% rewarding shape equally presented on the left and the right side (max 5
% diff)
pos_rewarding = zeros(nblck,ntrl);
for iblck = 1:nblck
    pos_rewarding(iblck,expe.rslt(iblck).feedback==1) = expe.rslt(iblck).respkb(expe.rslt(iblck).feedback==1);
    pos_rewarding(iblck,expe.rslt(iblck).feedback~=1) = 3-expe.rslt(iblck).respkb(expe.rslt(iblck).feedback~=1);
end
h = hist(pos_rewarding',3);
h = h(1:2,:);
dh = diff(h);
passed_eq_rew = all(abs(dh(1:2))<=5);

end