% FormatTrans transforms the 3D profile to 2D states and meas.
% Suppose all targets have same lifetime, or the number of steps 
% 'nStep' is constant for all targets

function [state, meas, nStep] = ...
	FormatTrans(profile, nTarg, densClt, Pd, r, noPlot)

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
poissClt = densClt*(xMax-xMin)*(yMax-yMin);
for i = 1 : nStep
	thisMeas = [];
	if i == 1
		% meas{1} is accurate (for initialization)
		for t = 1 : nTarg
			thisMeas = [thisMeas, state{t}([1 3], i)];
		end
	else
		for t = 1 : nTarg
			if rand < Pd % detected
				thisMeas = [thisMeas, ...
					state{t}([1 3], i) + normrnd(0, r, 2, 1)];
			end
		end

		nClt = poissrnd(poissClt);
		% make sure nClt > 0 to avoid the case that meas{i} is empty
		while nClt == 0
			nClt = poissrnd(poissClt);
		end
		cltMeas = [unifrnd(0.9*xMin, 1.1*xMax, 1, nClt); ...
			unifrnd(0.9*yMin, 1.1*yMax, 1, nClt)];
		thisMeas = [thisMeas, cltMeas];
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
		plot(state{t}(1, :), state{t}(3, :),'ro');
		hold on
	end
	for i = 1 : nStep
		measMat = meas{i};
		if size(measMat, 2) > 0
			plot(measMat(1, :), measMat(2, :), '*');
			hold on
		end
	end
end

if noPlot == 0
    hold off
end

