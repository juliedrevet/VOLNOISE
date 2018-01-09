function combisubj = genjitters(subj,nrev)
% this function generate squares of the 9 possible iti-stout jitter 
% combinations for subject subj for nrev number of reversals

% we want iti and stout jitters to be always different the 3 trials before 
% and 3 trials after a reversal for the 3 blocks so we use a dim 3 latin square
lsq = latin_square; % 12 possible latin squares for iti, 12 for stout

%% generate all possible iti-stout combination squares (12*12 = 144 squares)
allcombi = zeros(3,3,144); 
% logical array containing 1 if iti+stout combination are complementary 
% (i.e. if the square contains all 9 combinations)
tokeep   = false(1,144);
l = 1;
for i = 1:12
    for j = 1:12
        allcombi(:,:,l) = give_combination(lsq(:,:,i) ,lsq(:,:,j));
        tokeep(l)       = iscomplementary(lsq(:,:,i) ,lsq(:,:,j));
        l = l+1;
    end
end

%% only keep combination squares containing all 1:9 combinations
allcombi = allcombi(:,:,tokeep);
% reshape all possible iti-stout combination squares as for each 12 stout
% square, only 6 iti are complementary
allcombibis = zeros(3,3,6,12);
for k = 1:12
    allcombibis(:,:,:,k) = allcombi(:,:,((k-1)*6+1):k*6);
end

%% assign the first stout according to subject number
first_idx_stout = mod(subj-1,12)+1;

idx_stout = randperm(12,nrev); % choose a different stout latin square for each reversal
while idx_stout(1)~=first_idx_stout
    idx_stout = randperm(12,nrev);
end
idx_iti = randperm(6,1); % find complementary iti

combisubj(:,:,1) = allcombibis(:,:,idx_iti,idx_stout(1));

%% generate iti-stout combination squares for each following reversal
nperm_stout = 5e3;
nperm_iti   = 5e2;

found = false;
for q = 1:nperm_stout
    idx_stout = randperm(12,nrev);
    while idx_stout(1)~=first_idx_stout
        idx_stout = randperm(12,nrev);
    end
    for p = 1:nperm_iti
        for k = 2:nrev % for each reversal
            idx_iti = randperm(6,1); % find complementary iti
            combisubj(:,:,k) = allcombibis(:,:,idx_iti,idx_stout(k));
            % first column and 3 column should contain 1:9 combinations
            % different from previous square
            while sum(ismember(combisubj(:,1,k),combisubj(:,1,k-1)))==3 ||...
                    sum(ismember(combisubj(:,3,k),combisubj(:,3,k-1)))==3
                idx_iti = randperm(6,1);
                combisubj(:,:,k) = allcombibis(:,:,idx_iti,idx_stout(k));
            end        
        end
        % if all 9 combinations cannot be present in the 1st or last column of the latin square 
        % as there are too few reversals they must at least be all different
        if (nrev*3)<9 
            if (sum(ismember(1:9,reshape(combisubj(:,1,:),1,3*nrev)))>=(nrev*3))&&...
                    (sum(ismember(1:9,reshape(combisubj(:,3,:),1,3*nrev)))>=(nrev*3))
                found = true;
                break;
            end
        else % if all 9 combinations can be present, they should be equally distributed in 1st|last column
            h1 = hist(reshape(combisubj(:,1,:),1,3*nrev),1:9); % first column
            h2 = hist(reshape(combisubj(:,3,:),1,3*nrev),1:9); % last column
            if ~any(h1<floor(nrev/3))&&~any(h1>ceil(nrev/3))&&...
               ~any(h2<floor(nrev/3))&&~any(h2>ceil(nrev/3))
                found = true;
                break;
            end
        end
        if found
            break;
        end
    end
    if found
        break;
    end
end
if ~found
    error('unable to generate jitters')
end


end

function itistout = give_combination(iti,stout)
% return iti-stout jitters combinations (1:9) given iti and stout arrays
% (size not necessarily 3*3)

% 9 possible combinations
%          ITI STOUT
combi = [11; ... %1
         12; ... %2
         13; ... %3
         21; ... %4
         22; ... %5
         23; ... %6
         31; ... %7
         32; ... %8
         33];    %9
c = cell(size(iti));
itistout = zeros(size(iti));

for j = 1:size(iti,1)
    for k = 1:size(iti,2)
        c{j,k}=sprintf('%d%d',iti(j,k),stout(j,k));
        itistout(j,k) = find(str2double(c{j,k})==combi);
    end
end

end

function compl = iscomplementary(l1,l2)
% compare 2 latin squares and return 1 if they are complementary 
% (i.e. all 9 combinations fullfullied)

c = give_combination(l1,l2);

if length(unique(c)) == 9
    compl = true;
else
    compl = false;
end
end


