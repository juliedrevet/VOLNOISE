function [expe] = gen_prac(nprac)
%  GEN_PRAC  Generate VOLNOISE training experiment
%
%  Usage: [expe] = GEN_PRAC(nprac)
%
%  nprac: number of practice blocks

addpath ./Toolboxes/Rand
addpath ./Toolboxes/IO
addpath ./Toolboxes/Stimuli/

% check input arguments
if nargin < 1
    nprac = 1;
end

for i = 1:nprac
    fprintf('Generating practice block %d/%d...\n',i,nprac);
    blck = gen_blck(0);
    expe.blck(i) = blck;
end

end