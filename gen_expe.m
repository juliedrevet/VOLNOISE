function [expe] = gen_expe(subj,session)
%  GEN_EXPE  Generate VOLNOISE experiment
%
%  Usage: [expe] = GEN_EXPE(subj,session)
%
%  subj    - subject number
%  session - 1 or 2
%  volatility permutations are counter-balanced after group of 12 subjects.

addpath ./Toolboxes/Rand
addpath ./Toolboxes/IO
addpath ./Toolboxes/Stimuli/

% check input arguments
if nargin < 1
    error('Missing subject number!');
elseif nargin == 1
    error('Missing session!');
elseif ~isscalar(subj) || mod(subj,1) ~= 0
    error('Invalid subject number!');
elseif session ~= 1 && session ~= 2
    error('Invalid session number!'); 
end

%% experiment configuration parameters:
% volatility levels permutation
runs = latin_square; % dimension 3 (3 volatilities per block)
if session == 1
    volperm = runs(:,:,mod(subj,12)+1); % line: block, column: volatility level
else % session 2
    volperm = runs(:,:,mod(subj+4,12)+1); % all volatility level different
end


%% generate experiment blocks
for i = 1:size(volperm,2)
    fprintf('Generating block %d/%d...\n',i,size(volperm,2));
    blck = gen_blck(volperm(i,:));
    expe.blck(i) = blck;
end

end