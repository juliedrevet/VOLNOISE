function compl = iscompl(l1,l2)
% compare 2 latin squares and return 1 if they are complementary (all 9
% combinations fullfullied)
% verifiy complementarity
c = {};
m = 1;
for j = 1:3
    for k = 1:3
        c{m}=sprintf('%d%d',l1(j,k),l2(j,k));
        m = m+1;
    end
end
if length(unique(c)) == 9
    compl = true;
else
    compl = false;
end


end