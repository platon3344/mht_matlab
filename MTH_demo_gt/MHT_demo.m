% MHT_demo
% 代码中merty算法产生最优的假设，也就是新息值最小的前（）个假设
% prune剪枝算法，执行航迹产生和消亡。
clc;clear;close all;
noPlot = 0;
load CurveOne.mat profile;
tic;
% values of nTarg and T are determined by .mat file
nTarg = 1;
T = 2;
% assign values of other necessary parameters
densClt = 2e-7;
Pd = 0.9;
q = 500;
r = 50;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FormatTrans transforms the 3D profile to 2D states and meas.
% Suppose all targets have same lifetime, or the number of steps
% 'nStep' is constant for all targets
state = cell(1, nTarg);
xMax = -Inf; xMin = Inf;
yMax = -Inf; yMin = Inf;
for i = 1 : nTarg
    targ = profile{i};
    stateMat = (targ(:, [1 4 2 5]))';%x,vx,y,vy;
    state{i} = stateMat;
    xMax = max(xMax, max(stateMat(1, :)));
    xMin = min(xMin, min(stateMat(1, :)));
    yMax = max(yMax, max(stateMat(3, :)));
    yMin = min(yMin, min(stateMat(3, :)));
end

nStep = size(stateMat, 2);
meas = cell(1, nStep);
poissClt = densClt*(xMax-xMin)*(yMax-yMin);%poissClt为泊松系数，根据泊松杂波密度，产生泊松系数lamda
for i = 1 : nStep
    thisMeas = [];
    if i == 1
        % meas{1} is accurate (for initialization)
        for t = 1 : nTarg
            thisMeas = [thisMeas, state{t}([1 3], i)];%第一帧数据取真实值作为其实值，取第一行第三行，第i列数据
        end
    else
        for t = 1 : nTarg
            if rand < Pd % detected
                thisMeas = [thisMeas, state{t}([1 3], i) + normrnd(0, r, 2, 1)];%第二帧数据，开始根据随机产生添加噪声的测量数据
            end
        end
        
        nClt = poissrnd(poissClt);%根据泊松系数，产生泊松杂波个数
        % make sure nClt > 0 to avoid the case that meas{i} is empty
        while nClt == 0
            nClt = poissrnd(poissClt);
        end
        cltMeas = [unifrnd(0.9*xMin, 1.1*xMax, 1, nClt); unifrnd(0.9*yMin, 1.1*yMax, 1, nClt)];%x,y方向上产生1行nClt列的均匀分布矩阵
        thisMeas = [thisMeas, cltMeas];%将带噪声的数据和杂波，联合产生数据作为最终测量数据
    end
    meas{i} = thisMeas;
end

% the philosophy of this step is that N points seperate N-1
% parts. The first meas is for initialization, not estimation.
nStep = nStep - 1;

% plot if allowed
if noPlot == 0
    figure;
    for t = 1 : nTarg
        plot(state{t}(1, :), state{t}(3, :),'ro');%真实的测量值
        hold on
    end
    for i = 1 : nStep
        measMat = meas{i};
        if size(measMat, 2) > 0%每一列代表杂波或目标，行代表测量值x,y，大于0表示存在目标或者杂波
            plot(measMat(1, :), measMat(2, :), '*');
            hold on
        end
    end
end

if noPlot == 0
    hold off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = 3; % number of most probable hypos listed in each scan，有多少颗树，也就是有多个航迹族
N = 3; % height of the hypo tree. It seems scan depth = N-1，每颗树的最大深度
densNew = 0;

numTarget = nTarg;
%%initialization
F = [1 T 0 0; 0 1 0 0; 0 0 1 T; 0 0 0 1];
H = [1 0 0 0; 0 0 1 0];
P = 0; % zero initial covariance
Q = q*[T^3/3 T^2/2 0 0; T^2/2 T 0 0; 0 0 T^3/3 T^2/2; 0 0 T^2/2 T];
R = diag([r^2, r^2]);

maxNumOfHypo = 0; % max number of hypoes kept at the same time
for i = 1 : N
    maxNumOfHypo = maxNumOfHypo + M^i;%不知道该地方如何计算，为什么是3^1 + 3^2 + 3^3 = 39 ，暂时先放一放 
end

% fill in 1st cell in cellHypo: {{asso prob}  }
% note: the value of root prob is not important
cellHypo = {{(1:numTarget)' 0}};%产生2个cell，分别赋值第一个cell为目标序列数组，第二个cell为0

% fill in 1st cell in cellTarg: {{numTarg*{idx lifePoint X P}}  }
maxLifePoint = 3;
cellTmp = cell(1, numTarget); % number of target is known and constant
for i = 1 : numTarget
    tmp = state{i};
    cellTmp{i} = {i maxLifePoint tmp(:, 1) P};%这个为每个目标的数据结构
end
cellTarg = {cellTmp};%目标cell组，因为只有一个目标所以1*1

% fill in 1st cell in cellEstm: {numTarg*{idx startTime matX}}
cellEstm = cell(1, numTarget);
for i = 1 : numTarget
    tmp = state{i};
    cellEstm{i} = {i 0 tmp(:, 1)};%赋初值
end

%%evolution stage
head = 1;
rear = 1;

if noPlot == 0
    figure;
end

for t = 1 : nStep
    % hypoes and estims formed last step as seeds to generate new ones
    
    cellHypoSeed = {cellHypo{head:rear}};
    cellTargSeed = {cellTarg{head:rear}};
    
    % predict stage of Kalman filter
    cellTargSeed = KF_MHT_Predict(cellTargSeed, F, Q);%通过预测方程，得到预测值位置更新，和预测协方差矩阵
    
    % generate new hypoes
    Z = meas{t+1}; % get measurement of this time step%测到测量值坐标，包含了目标和杂波
    if isempty(Z)
        error('There in no measurement at step %d.', t);
    else
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % GenHypo lists M most possible hypotheses, calculates their
        %假设产生操作
        len = size(cellHypoSeed, 2);%获取当前假设的个数
        cellHypoNew = [];%新的假设
        for i = 1 : len
            targTemp = cellTargSeed{i}; % {nTarg * {idx, lifePoint, X, P}}
            % in accord with cellHypoSeed{i}
            maxTargIdx = max(cellfun(@(v) v{1}, targTemp));%目标索引号最大值
            nMeas = size(Z, 2);%目标和杂波总个数
            
            % form the log-prob
            % matrix，计算量测值和每个航迹之间的关联概率矩阵，作为代价矩阵，用匈牙利算法和murty算法
            ProbMat = GenProbMat(maxTargIdx, targTemp,Z, H, R, Pd, densNew, densClt);
            
            % use Murty's algorithm to find out M minimum assignments
            cellHypoTmp = Murty(ProbMat, M);
            
            hypoSeed = cellHypoSeed{i};
            p1 = hypoSeed{2}; % log-prob of this hypoSeed
            nT = sum(hypoSeed{1} > 0); % number of targets deteced in hypoSeed
            p2 = -nT * log(1-Pd);
            for j = 1 : M
                hypo = cellHypoTmp{j};%遍历假设
                asso = hypo{1};
                
                % set indicators of clutter to 0 in assignment matrix
                % see mat1, mat2 and mat3 in GenProbMat
                asso(asso > maxTargIdx + nMeas) = 0;%代码是否应该是asso(asso > maxTargIdx) = 0
                
                % make indicators to new targets more compact
                idxNew = find(asso > maxTargIdx);
                asso(idxNew) = maxTargIdx+1:maxTargIdx+length(idxNew);
                
                % go on with probability calculation
                prob = hypo{2};
                prob = prob + p1 + p2;
                
                hypo{1} = asso;
                hypo{2} = prob;
                cellHypoTmp{j} = hypo;
            end
            cellHypoNew = [cellHypoNew, cellHypoTmp];%cellHypoNew得到航迹和量测值的最终关联矩阵，然后利用卡尔曼更新
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % update stage of Kalman filter
        % 根据假设更新每一个目标的位置
        cellTargNew = KF_MHT_Update(cellTargSeed, cellHypoNew,M, Z, H, R, maxLifePoint);%卡尔曼更新
        
        % expand cellTarg and cellHypo. NOTE THE DIFFERENCE BETWEEN {} & []
        cellHypo = [cellHypo, cellHypoNew];
        cellTarg = [cellTarg, cellTargNew];
        
        % prune and update剪枝，最大深度为N
        if t < N%1到2帧不需要剪枝，因为都没有形成深度为3航迹
            head = rear + 1;
            rear = rear + M^t;
        else % positions of head and rear are now fixed
            [cellEstm, cellHypo, cellTarg] = Prune(cellEstm, cellHypo, cellTarg, M, N, t, maxLifePoint);%剪枝叶操作
        end
    end
    
    for i = 1 : size(cellEstm, 2)
        stateEstm = cellEstm{i}{3};
        if noPlot == 0
            plot(state{i}(1, 1:t+1), state{i}(3, 1:t+1), '-', stateEstm(1, :), stateEstm(3, :), '*');
            hold on;
            drawnow;
        end
    end
end

if noPlot == 0
    hold off;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1 : length(cellEstm)
    cellEstm{i} = cellEstm{i}{3}; % extract the state
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate mean error
first = 2; % estmMHT includes accurate initial state
last = (nStep+1) - (N-1);
[errRMS, lose] = Analyse(first, last, cellEstm, state, nTarg);
rstMHT = [errRMS lose]; % MHT tracking results
disp('    err_x     err_vx     err_y     err_vy     lose  sum');
disp(rstMHT);
disp(sum(errRMS));
fprintf('MHT DONE\n\n');
toc;