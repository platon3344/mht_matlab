% KF_MHT_PREDICT  Perform Kalman Filter prediction step in MHT algorithm
%
% Syntax:
%   [X,P] = KF_PREDICT(X,P,A,Q,B,U)
%
% In:
%   cellTargSeed - a cell array containing description of targets. 
%       Each cell contains information of a full collection of targets
%       derived according to the corresponding hypo. One target was
%       described with {idx, lifePoint, X, P}. 
%   A - Transition matrix of discrete model (optional, default identity)
%   Q - Process noise of discrete model     (optional, default zero)
%   B - Input effect matrix                 (optional, default identity)
%   U - Constant input                      (optional, default empty)
%
% Out:
%  cellTargSeed - with X, P of each target in each cell modified
%   
% Description:
%   Perform Kalman Filter prediction step. The model is
%
%     x[k] = A*x[k-1] + B*u[k-1] + q,  q ~ N(0,Q).
% 
%   The predicted state is distributed as follows:
%   
%     p(x[k] | x[k-1]) = N(x[k] | A*x[k-1] + B*u[k-1], Q[k-1])
%
%   The predicted mean x-[k] and covariance P-[k] are calculated
%   with the following equations:
%
%     m-[k] = A*x[k-1] + B*u[k-1]
%     P-[k] = A*P[k-1]*A' + Q.
%
%   If there is no input u present then the first equation reduces to
%     m-[k] = A*x[k-1]

function cellTargSeed = KF_MHT_Predict(cellTargSeed, A, Q, B, U)

if nargin < 2
    A = [];
end
if nargin < 3
    Q = [];
end
if nargin < 4
    B = [];
end
if nargin < 5
    U = [];
end

for i = 1 : size(cellTargSeed, 2)
    oneCell = cellTargSeed{i};
    for j = 1 : size(oneCell, 2)
        % if lifePoint == 0, then pass without processing
        if oneCell{j}{2} == 0
            continue;
        end

        X = oneCell{j}{3};
        P = oneCell{j}{4};

        % Apply defaults
        if isempty(A)
            A = eye(size(X,1));
        end
        if isempty(Q)
            Q = zeros(size(X,1));
        end
        if isempty(B) && ~isempty(U)
            B = eye(size(X,1),size(U,1));
        end

        % Perform prediction
        if isempty(U)
            X = A * X;
            P = A * P * A' + Q;
        else
            X = A * X + B * U;
            P = A * P * A' + Q;
        end

        oneCell{j}{3} = X;
        oneCell{j}{4} = P;
    end
    cellTargSeed{i} = oneCell;
end
