function itistout = give_combi(iti,stout)
% 9 possible combinations
%          ITI STOUT
combi = [11; ...
    12; ...
    13; ...
    21; ...
    22; ...
    23; ...
    31; ...
    32; ...
    33];
c1 = {};
c2 = {};
itistout = zeros(size(iti));

for j = 1:size(iti,1)
    for k = 1:size(iti,2)
        c1{j,k}=sprintf('%d%d',iti(j,k),stout(j,k));
        itistout(j,k) = find(str2num(c1{j,k})==combi);
        c2{j,k}=sprintf('%d%d',iti(j,k),stout(j,k));
        itistout(j,k) = find(str2num(c2{j,k})==combi);
    end
end

end