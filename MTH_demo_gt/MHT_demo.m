% MHT_demo
% ������merty�㷨�������ŵļ��裬Ҳ������Ϣֵ��С��ǰ����������
% prune��֦�㷨��ִ�к���������������
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
poissClt = densClt*(xMax-xMin)*(yMax-yMin);%poissCltΪ����ϵ�������ݲ����Ӳ��ܶȣ���������ϵ��lamda
for i = 1 : nStep
    thisMeas = [];
    if i == 1
        % meas{1} is accurate (for initialization)
        for t = 1 : nTarg
            thisMeas = [thisMeas, state{t}([1 3], i)];%��һ֡����ȡ��ʵֵ��Ϊ��ʵֵ��ȡ��һ�е����У���i������
        end
    else
        for t = 1 : nTarg
            if rand < Pd % detected
                thisMeas = [thisMeas, state{t}([1 3], i) + normrnd(0, r, 2, 1)];%�ڶ�֡���ݣ���ʼ�������������������Ĳ�������
            end
        end
        
        nClt = poissrnd(poissClt);%���ݲ���ϵ�������������Ӳ�����
        % make sure nClt > 0 to avoid the case that meas{i} is empty
        while nClt == 0
            nClt = poissrnd(poissClt);
        end
        cltMeas = [unifrnd(0.9*xMin, 1.1*xMax, 1, nClt); unifrnd(0.9*yMin, 1.1*yMax, 1, nClt)];%x,y�����ϲ���1��nClt�еľ��ȷֲ�����
        thisMeas = [thisMeas, cltMeas];%�������������ݺ��Ӳ������ϲ���������Ϊ���ղ�������
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
        plot(state{t}(1, :), state{t}(3, :),'ro');%��ʵ�Ĳ���ֵ
        hold on
    end
    for i = 1 : nStep
        measMat = meas{i};
        if size(measMat, 2) > 0%ÿһ�д����Ӳ���Ŀ�꣬�д������ֵx,y������0��ʾ����Ŀ������Ӳ�
            plot(measMat(1, :), measMat(2, :), '*');
            hold on
        end
    end
end

if noPlot == 0
    hold off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = 3; % number of most probable hypos listed in each scan���ж��ٿ�����Ҳ�����ж��������
N = 3; % height of the hypo tree. It seems scan depth = N-1��ÿ������������
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
    maxNumOfHypo = maxNumOfHypo + M^i;%��֪���õط���μ��㣬Ϊʲô��3^1 + 3^2 + 3^3 = 39 ����ʱ�ȷ�һ�� 
end

% fill in 1st cell in cellHypo: {{asso prob}  }
% note: the value of root prob is not important
cellHypo = {{(1:numTarget)' 0}};%����2��cell���ֱ�ֵ��һ��cellΪĿ���������飬�ڶ���cellΪ0

% fill in 1st cell in cellTarg: {{numTarg*{idx lifePoint X P}}  }
maxLifePoint = 3;
cellTmp = cell(1, numTarget); % number of target is known and constant
for i = 1 : numTarget
    tmp = state{i};
    cellTmp{i} = {i maxLifePoint tmp(:, 1) P};%���Ϊÿ��Ŀ������ݽṹ
end
cellTarg = {cellTmp};%Ŀ��cell�飬��Ϊֻ��һ��Ŀ������1*1

% fill in 1st cell in cellEstm: {numTarg*{idx startTime matX}}
cellEstm = cell(1, numTarget);
for i = 1 : numTarget
    tmp = state{i};
    cellEstm{i} = {i 0 tmp(:, 1)};%����ֵ
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
    cellTargSeed = KF_MHT_Predict(cellTargSeed, F, Q);%ͨ��Ԥ�ⷽ�̣��õ�Ԥ��ֵλ�ø��£���Ԥ��Э�������
    
    % generate new hypoes
    Z = meas{t+1}; % get measurement of this time step%�⵽����ֵ���꣬������Ŀ����Ӳ�
    if isempty(Z)
        error('There in no measurement at step %d.', t);
    else
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % GenHypo lists M most possible hypotheses, calculates their
        %�����������
        len = size(cellHypoSeed, 2);%��ȡ��ǰ����ĸ���
        cellHypoNew = [];%�µļ���
        for i = 1 : len
            targTemp = cellTargSeed{i}; % {nTarg * {idx, lifePoint, X, P}}
            % in accord with cellHypoSeed{i}
            maxTargIdx = max(cellfun(@(v) v{1}, targTemp));%Ŀ�����������ֵ
            nMeas = size(Z, 2);%Ŀ����Ӳ��ܸ���
            
            % form the log-prob
            % matrix����������ֵ��ÿ������֮��Ĺ������ʾ�����Ϊ���۾������������㷨��murty�㷨
            ProbMat = GenProbMat(maxTargIdx, targTemp,Z, H, R, Pd, densNew, densClt);
            
            % use Murty's algorithm to find out M minimum assignments
            cellHypoTmp = Murty(ProbMat, M);
            
            hypoSeed = cellHypoSeed{i};
            p1 = hypoSeed{2}; % log-prob of this hypoSeed
            nT = sum(hypoSeed{1} > 0); % number of targets deteced in hypoSeed
            p2 = -nT * log(1-Pd);
            for j = 1 : M
                hypo = cellHypoTmp{j};%��������
                asso = hypo{1};
                
                % set indicators of clutter to 0 in assignment matrix
                % see mat1, mat2 and mat3 in GenProbMat
                asso(asso > maxTargIdx + nMeas) = 0;%�����Ƿ�Ӧ����asso(asso > maxTargIdx) = 0
                
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
            cellHypoNew = [cellHypoNew, cellHypoTmp];%cellHypoNew�õ�����������ֵ�����չ�������Ȼ�����ÿ���������
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % update stage of Kalman filter
        % ���ݼ������ÿһ��Ŀ���λ��
        cellTargNew = KF_MHT_Update(cellTargSeed, cellHypoNew,M, Z, H, R, maxLifePoint);%����������
        
        % expand cellTarg and cellHypo. NOTE THE DIFFERENCE BETWEEN {} & []
        cellHypo = [cellHypo, cellHypoNew];
        cellTarg = [cellTarg, cellTargNew];
        
        % prune and update��֦��������ΪN
        if t < N%1��2֡����Ҫ��֦����Ϊ��û���γ����Ϊ3����
            head = rear + 1;
            rear = rear + M^t;
        else % positions of head and rear are now fixed
            [cellEstm, cellHypo, cellTarg] = Prune(cellEstm, cellHypo, cellTarg, M, N, t, maxLifePoint);%��֦Ҷ����
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