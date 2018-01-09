function compl = find_complementary(lsq)

for i = 1:12
    l = 1;
    for j = 1:12
        if iscompl(lsq(:,:,mod(i-1,12)+1),lsq(:,:,mod(j-1,12)+1))
            compl(i,l) = j;
            l = l+1;
        end
    end
end
