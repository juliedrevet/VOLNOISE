function [passed,match_outcome] = response_feedback(expe)
% compares the subject responses and the received feedback to contruct back 
% the rewarding sequence (only possible when subject response exists!). 
% This sequence is then compared to the one that should have been presented (expe.blck.outcome)

ntrl  = length(expe.rslt(1).resp);
nblck = length(expe.rslt);
constr_outcome = zeros(nblck,ntrl);
match_outcome  = zeros(nblck,ntrl); % 1 if constructed outcome matches blck.outcome, 0 otherwise
noresp = 0;

for iblck = 1:nblck
    constr_outcome(iblck,1:ntrl) = expe.rslt(iblck).resp;
    % opposite answer if feedback was 0
    constr_outcome(iblck,(expe.rslt(iblck).feedback(1:ntrl)==0)) = 3-constr_outcome(iblck,(expe.rslt(iblck).feedback(1:ntrl)==0));
    match_outcome(iblck,1:ntrl) = constr_outcome(iblck,1:ntrl)==expe.blck(iblck).outcome(1:ntrl);
    % count the number of times the subject did not respond
    noresp = noresp+sum(expe.rslt(iblck).resp(1:ntrl) == 0);
end

% verification test passed if the reconstructed outcome matches the
% blck.outcome structure (except for the trials the subject did not answer
% as nothing can be inferred about the correctness of the response)
passed = (sum(match_outcome(:))==(nblck*ntrl)-noresp);
end