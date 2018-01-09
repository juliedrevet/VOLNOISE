%% verification of experiment presentation

% load all experiments (TODO: all subj? one subj?)
[FileName,PathName] = uigetfile('*.mat','Select subject experiment structure');
expe = importdata(fullfile(PathName,FileName));

%
fprintf('Reconstructed rewarding sequence matches generated sequence...')
[passed_outcome,match_outcome] = response_feedback(expe);
if passed_outcome
    fprintf('OK\n')
else
    fprintf('NO!\n')
end
