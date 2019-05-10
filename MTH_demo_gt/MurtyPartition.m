% MurtyPartition partitioin node N with its minimum assignment a
% 
% Syntax: 
%   nodeList = MurtyPartition(N, a, type)
%   
% In: 
%   N - in Murty's original paper, N is a "node", i.e. a non empty 
%       subset of A, which contains all assignment schemes. Here N is
%       expressed with a cell {Inclu, Exclu}, see nodeList below.
%   a - a nMeas*1 vector containing one assignment scheme.
%   type - type == 0 for N-to-N assignment problem, type == 1 for 
%       M-to-N assignment problem, where M > N, e.g. assign M jobs to 
%       N worker. 
%   
% Out: 
%   nodeList - a cell array containing the list of partition of N. The
%       union of all assignments to all partitions and assignment {a} 
%       forms a complete set of assignments to N. One partition is ex-
%       pressed with a cell {Inclu, Exclu}, where
%         Inclu is a ?*2 matrix containing entries must be included in 
%           the assignments for this partition. Each row [idx, assign];
%         Exclu is a ?*2 matrix containing entries must be excluded in 
%           the assignments for this partition.

function nodeList = MurtyPartition(N, a, type)

nMeas = size(a, 1); 
a = [(1:nMeas)', a]; % add a column of index
a1 = intersect(N{1}, a, 'rows'); % entries must be included
a2 = setdiff(a, a1, 'rows'); % interchangeable entries

if type == 0 % N-to-N assignment
    nodeList = cell(1, size(a2, 1)-1);
else % M-to-N assignment, M > N
    nodeList = cell(1, size(a2, 1));
end

for i = 1 : length(nodeList)
    if i == 1
        Inclu = N{1};
    else
        Inclu = [N{1}; a2(1:i-1, :)];
    end
    Exclu = [N{2}; a2(i, :)];
    tmp = cell(1, 2);
    tmp{1} = Inclu;
    tmp{2} = Exclu;
    nodeList{i} = tmp;
end




