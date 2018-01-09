function opp = find_opposite(lsq)


for i = 1:12
    k = 1;
    for j = 1:12
        if ~any(lsq(:,:,i) == lsq(:,:,j))
            opp(i,k) = j;
            k = k+1;
        end
    end
end
end
