% KF_MHT_Update generates target states according to multiple hypotheses
% 
% Syntax:
%   cellTargNew = KF_MHT_Update(cellTargSeed, cellHypoNew,...
%       M, Z, H, R, maxLifePoint)
%   
% In: 
%   cellTargSeed - a cell array containing description of targets. 
%       Each cell contains information of a full collection of targets
%       derived according to the corresponding hypo. One target was
%       described with {idx, lifePoint, X, P}. 
%   cellHypoNew - a cell array containing newly generated hypoes and 
%       their probabilities. 
%       Note: size(cellHypoNew, 2) == M * size(cellTargSeed, 2)
%   M - number of hypoes generated using Murty's algorithm
%   Z - measurements at this time step
%   H - measurement matrix
%   R - measurement covariance
%   maxLifePoint - maximum life point: each target maintains a "life 
%       point", at each time step, if this target is not detected, then
%       subtract 1 from its life point; else if its life point is less
%       than maxLifePoint, then add 1 to its life point.
%   
% Out: 
%   cellTargNew - a cell array containing newly generated targets info
%       that is one-one correspondent to cellHypoNew. Each cell in 
%       cellHypoNew leads to one cell in cellTargNew which gives a full
%       set of state estimations for all the targets.  
%
% Description: 
%  cellHypo{i} is {asso_i, prob_i}
%  cellTarg{i} is {nTarg_i*{targInfo}}, 
%      and targInfo{i} is {idx_i, lifePoint_i, X_i, P_i}

function cellTargNew = KF_MHT_Update(cellTargSeed, cellHypoNew,...
    M, Z, H, R, maxLifePoint)

rear = 0;
cellTargNew = cell(1, size(cellHypoNew, 2));

for i = 1 : size(cellTargSeed, 2)
    head = rear + 1;
    rear = rear + M;
    for j = head : rear
        % use cellHypoNew{j} to update cellTargSeed{i}, thus form
        % cellTargNew{j}
        asso = cellHypoNew{j}{1};
        nExistedTarg = size(cellTargSeed{i}, 2); 
        maxTargIdx = max(cellfun(@(v) v{1}, cellTargSeed{i}));
        nNewTarg = sum(asso > maxTargIdx);
        
        aCase = cell(1, nExistedTarg+nNewTarg);
        % 1. deal with existed (in last step) targets 
        for k = 1 : nExistedTarg % for each target
            aTarg = cellTargSeed{i}{k}; % one target
            idx = aTarg{1}; % the index of aTarg
            lifePoint = aTarg{2}; 
            X = aTarg{3};
            P = aTarg{4};
            
            if lifePoint == 0 % a disappeared target
                aCase{k} = aTarg;
                continue; % just pass it
            end
            flg = find(asso == idx);
            if isempty(flg) % there is no meas asso with aTarg
                lifePoint = lifePoint - 1;
            else
                aMeas = Z(:, flg); % the meas asso with aTarg
                
                % Kalman filter update stage
                innov = aMeas - H*X;
                S = R + H*P*H';
                G = P*H'/S;
                X = X + G*innov;
                P = P - G*S*G';
                
                if lifePoint < maxLifePoint
                    lifePoint = lifePoint + 1;
                end
            end
                
            aTarg{2} = lifePoint;
            aTarg{3} = X;
            aTarg{4} = P;
           
            aCase{k} = aTarg;
        end
        % 2. deal with newly observed tentative targets
        for k = 1 : nNewTarg
            idx = maxTargIdx + k;
            flg = find(asso == idx);
            aMeas = Z(:, flg); 
            
            % initialize a new target
            aTarg = cell(1, 4);
            aTarg{1} = idx;
            aTarg{2} = 1; % lifePoint
            aTarg{3} = [aMeas(1), 0, aMeas(2), 0]'; % X
            aTarg{4} = diag([1 1 1 1]); % P
            
            aCase{idx} = aTarg;
        end
    cellTargNew{j} = aCase;
    end
end





