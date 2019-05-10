% Murty's algorithm finds out the kth minimum assignments, k = 1, 2, ...
% 
% Syntax: 
%   solution = Murty(costMat, k)
%   
% In: 
%   costMat - nMeas*nTarg cost matrix.
%   k - the command number controlling the output size.
%   
% Out: 
%   solution - cell array containing the minimum, 2nd minimum, ..., 
%       kth minimum assignments and their costs. Each solution{i} 
%       contains {assgmt, cost} where assgmt is an nMeas*1 matrix 
%       giving the ith minimum assignment; cost is the cost of this 
%       assignment.
      
function solution = Murty(costMat, k)

solution = cell(1, k);

t = 1;
[assgmt, cost] = Hungarian(costMat);%行为测量值个数，列为航迹个数

solution{1} = {assgmt, cost}; 

% xxxRec stands for 'record'
nodeRec = cell(1, 2);
nodeRec{1} = [1,1];
assgmtRec = assgmt;

nodeList = MurtyPartition(nodeRec, assgmtRec, 1);

while t < k
    tmp = []; % structure space for temporary (assgmt, cost) storage
    minCost = Inf;
    idxRec = -1;

    % try to find one node in the nodeList with the minimum cost 
    for i = 1 : size(nodeList, 2)
        node = nodeList{i};
        Inclu = node{1};
        Exclu = node{2};
        mat = costMat;
        for j = 1 : size(Inclu, 1) % restrict: assignments must be included
            best = mat(Inclu(j, 1), Inclu(j, 2));
            mat(Inclu(j, 1), :) = Inf;
            mat(Inclu(j, 1), Inclu(j, 2)) = best;
        end
        for j = 1 : size(Exclu, 1) % restrict: assignments must be excluded
            mat(Exclu(j, 1), Exclu(j, 2)) = Inf;
        end
        
        [assgmt, cost] = Hungarian(mat);
        
        if ismember(0, assgmt) % cost have to be inf
            continue;
        elseif cost < minCost%寻找最小的代价函数值，对应的分配方法
            minCost = cost;
            nodeRec = node;
            assgmtRec = assgmt;
            idxRec = i;
        end
    end
    
    
    if idxRec == -1 % all node in the nodeList leads to inf cost
        for i = t+1 : k
            solution{i} = solution{t};
        end
        t = k;
    else    
        t = t + 1;
        solution{t} = {assgmtRec, minCost};
        idx = setdiff(1:size(nodeList, 2), idxRec);
        nodeList = [nodeList(idx), MurtyPartition(nodeRec, assgmtRec, 1)];
    end
end

