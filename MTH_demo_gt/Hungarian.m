% Hungarian deals with distmMat with negative elements

function [assgmt, cost] = Hungarian(distMat)

minEmt = min(distMat(:)); % minimum element in distMat

if minEmt >= 0
    [assgmt, cost] = assignmentoptimal(distMat);
else
    distMat = distMat - minEmt;
    [assgmt, cost] = assignmentoptimal(distMat);
    cost = cost + minEmt*size(distMat, 1);
end

% if ismember(0, assgmt)
%     error('All possible assignments leads to infinite cost.');
% end