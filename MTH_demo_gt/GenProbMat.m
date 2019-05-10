% GenProbMat forms the cost matrix, or log-prob matrix here
% 
% Syntax:
%     ProbMat = GenProbMat(maxTargIdx, targ, ...
%         Z, H, R, pd, densNew, densClt)
%     
% In: 
%   maxTargIdx - the maximum index of target that had been associtated 
%       with a measurement
%   targ - a cell array {nTarg * {idx, lifePoint, X, P}}. Each cell in 
%       targ is a cell describing one specific target.
%   Z - a dMeas*nMeas matrix containing meas on this step
%   H - measurement matrix
%   R - measurement covariance matrix
%   pd - probability of detection
%   densNew - density of new target in the surveillance space
%   densClt - density of clutter in the surveillance space  
%   
% Out: 
%   ProbMat - a  matrix containing the log probabilities of each 
%       assocition. 
      
function ProbMat = GenProbMat(maxTargIdx, targ, Z, H, R, pd, densNew, densClt)

nMeas = size(Z, 2);

% 1st part: meas associated with existed target 测量值和已经存在航迹(也就这个里面target)关联
mat1 = zeros(nMeas, maxTargIdx);
for i = 1 : nMeas
    for j = 1 : maxTargIdx
        % find the targ cell whose index is j
        thisTarg = targ{cellfun(@(v) v{1}, targ) == j};%循环targ的每一个cell，并取该cell的第一个函数的值，并判断该值是否等于j
        
        innov = H*thisTarg{3} - Z(:, i); % meas innovation计算新息
        S = H*thisTarg{4} * H' + R; % prior cov of innov计算协方差矩阵
        dim = length(thisTarg{3});%得到目标状态维度，x,vx,y,vy，一共四维
        
        if thisTarg{2} == 0 % lifePoint == 0, a disappeared targ
            mat1(i, j) = Inf;
        else
            x1 = log((1-pd)/pd); 
            x2 = 0.5*(innov'*inv(S)*innov + dim*log(2*pi) + log(det(S)));
            mat1(i, j) = x1 + x2;%计算航迹(targ)和每一个量测值的概率值
        end
    end
end

% 2nd part: meas associated with new target
mat2 = inf(nMeas, nMeas); 
for i = 1 : nMeas
    mat2(i, i) = -log(densNew);%计算每一个新的量测值产生概率
end

% 3rd part: meas associated with clutter
mat3 = inf(nMeas, nMeas);
for i = 1 : nMeas
    mat3(i, i) = -log(densClt);
end

% put them together
ProbMat = [mat1, mat2, mat3];