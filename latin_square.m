function squares = latin_square
% only for dim 3 (12 possibilities)

c1 = [1 2 3; ...
      2 3 1; ...
      3 1 2];
c2 = [1 3 2; ...
      2 1 3; ...
      3 2 1];
 
all_perms = perms(1:3);
squares   = zeros(3,3,12,'int8');

% first line
for k = 1:size(all_perms,1)
    squares(:,:,k) = c1(all_perms(k,:),:);
    squares(:,:,k+size(all_perms,1)) = c2(all_perms(k,:),:);
end

% re-organize position to facilitate subject attribution
% (session 1 ~ mod(12) || session 2 ~ mod(12) + 4 >> complete different square
sort_idx =  [1 10 12 3 5 11 9 2 4 7 8 6];
squares = squares(:,:,sort_idx);

% verify they are all unique
for k = 1:size(all_perms,1)
    s = squares(:,:,k);
    for j = 1:size(all_perms,1)
        if (s == squares(:,:,j))
            if k~=j
                sprintf('Square %d identical to square %d',k,j)
            end
        end
    end

end
end