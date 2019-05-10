function cellEstm = MHT(state, meas, T, nStep, numTarget, ...
    noPlot, densClt, densNew, Pd, q, r, M, N)

%% initialization
F = [1 T 0 0; 0 1 0 0; 0 0 1 T; 0 0 0 1];
H = [1 0 0 0; 0 0 1 0];
P = 0; % zero initial covariance
Q = q*[T^3/3 T^2/2 0 0; T^2/2 T 0 0; 0 0 T^3/3 T^2/2; 0 0 T^2/2 T];
R = diag([r^2, r^2]);

maxNumOfHypo = 0; % max number of hypoes kept at the same time
for i = 1 : N
    maxNumOfHypo = maxNumOfHypo + M^i;
end

% fill in 1st cell in cellHypo: {{asso prob}  }
% note: the value of root prob is not important
cellHypo = {{(1:numTarget)' 0}}; 

% fill in 1st cell in cellTarg: {{numTarg*{idx lifePoint X P}}  }
maxLifePoint = 3;
cellTmp = cell(1, numTarget); % number of target is known and constant
for i = 1 : numTarget
    tmp = state{i};
    cellTmp{i} = {i maxLifePoint tmp(:, 1) P};
end
cellTarg = {cellTmp};

% fill in 1st cell in cellEstm: {numTarg*{idx startTime matX}}
cellEstm = cell(1, numTarget); 
for i = 1 : numTarget
    tmp = state{i};
    cellEstm{i} = {i 0 tmp(:, 1)};
end

%% evolution stage
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
    cellTargSeed = KF_MHT_Predict(cellTargSeed, F, Q); 

    % generate new hypoes
    Z = meas{t+1}; % get measurement of this time step
	if isempty(Z)
		error('There in no measurement at step %d.', t);
	else
		cellHypoNew = GenHypo(cellHypoSeed, cellTargSeed, M,...
			Z, H, R, Pd, densNew, densClt);

		% update stage of Kalman filter
		cellTargNew = KF_MHT_Update(cellTargSeed, cellHypoNew,...
			M, Z, H, R, maxLifePoint);
    
		% expand cellTarg and cellHypo. NOTE THE DIFFERENCE BETWEEN {} & []
		cellHypo = [cellHypo, cellHypoNew];
		cellTarg = [cellTarg, cellTargNew];

		% prune and update
		if t < N
			head = rear + 1;
			rear = rear + M^t;
		else % positions of head and rear are now fixed
			[cellEstm, cellHypo, cellTarg] = Prune(cellEstm,...
				cellHypo, cellTarg, M, N, t, maxLifePoint);

		end
	end
	
	for i = 1 : size(cellEstm, 2)
		stateEstm = cellEstm{i}{3};

		if noPlot == 0
			plot(state{i}(1, 1:t+1), state{i}(3, 1:t+1), '-', ...
                stateEstm(1, :), stateEstm(3, :), 'o');
			hold on;
			drawnow;
		end
	end
end

if noPlot == 0
    hold off;
end
