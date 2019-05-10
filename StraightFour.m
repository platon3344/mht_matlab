function [state, meas, T, N, numTarget] = StraightFour(noPlot)

v1 = 0.009;
v2 = 0.009;
v3 = 0.009;
v4 = 0.0095;

T = 4;
N = 25; % simulation time: 25 * 4s
F = [1 T; 0 1];
q = 1e-17; % power spectral density of the process noise, Esti. P. 270
Q = q * [T^3/3 T^2/2; T^2/2 T];
r = 0.009; % measurement standard deviation
Pd = 0.9; % probability of detection
numTarget = 4;

clutterDens = 10; % clutter density
xMax = 1.1; xMin = -0.2;
yMax = 0.25; yMin = -0.2;
poissNum = clutterDens*(xMax-xMin)*(yMax-yMin);

state1 = [0 v1 0 0]'; % [x, vx, y, vy]
state2 = [0 v2 0.05 0]';
state3 = [0 v3 0.1 0]';
state4 = [0 v4*cos(deg2rad(20)) 0.2 -v4*sin(deg2rad(20))]';

meas = cell(N+1, 1);
state = cell(4, 1);
meas{1} = [state1(1) state2(1) state3(1) state4(1); ...
    state1(3) state2(3) state3(3) state4(3)];

for i = 1 : N
    stateX = F * state1(1:2, i) + DrawMultiNorm([0; 0], Q);
    stateY = F * state1(3:4, i) + DrawMultiNorm([0; 0], Q);
    state1 = [state1, [stateX; stateY]];
        
    stateX = F * state2(1:2, i) + DrawMultiNorm([0; 0], Q);
    stateY = F * state2(3:4, i) + DrawMultiNorm([0; 0], Q);
    state2 = [state2, [stateX; stateY]];

    stateX = F * state3(1:2, i) + DrawMultiNorm([0; 0], Q);
    stateY = F * state3(3:4, i) + DrawMultiNorm([0; 0], Q);
    state3 = [state3, [stateX; stateY]];
    
    stateX = F * state4(1:2, i) + DrawMultiNorm([0; 0], Q);
    stateY = F * state4(3:4, i) + DrawMultiNorm([0; 0], Q);
    state4 = [state4, [stateX; stateY]];    
    
    detect = find(unifrnd(0, 1, 4, 1) < Pd);
    posX = [state1(1,i+1) state2(1,i+1) state3(1,i+1) state4(1,i+1)];
    posY = [state1(3,i+1) state2(3,i+1) state3(3,i+1) state4(3,i+1)];
    measTarget = [posX(detect); posY(detect)]...
        + normrnd(0, r, 2, length(detect));
    numClutter = poissrnd(poissNum);
    measClutter = [unifrnd(xMin, xMax, 1, numClutter); ...
        unifrnd(yMin, yMax, 1, numClutter)];
    
    meas{i+1} = [measTarget, measClutter];
end

state{1} = state1;
state{2} = state2;
state{3} = state3;
state{4} = state4;

if noPlot == 0
    figure;
    plot(state1(1, :), state1(3, :), state2(1, :), state2(3, :), ...
         state3(1, :), state3(3, :), state4(1, :), state4(3, :));
    hold on
end

for i = 1 : N
    measMat = meas{i};
    if noPlot == 0 && size(measMat, 2) > 0
        plot(measMat(1, :), measMat(2, :), '*');
        hold on
    end
end

if noPlot == 0
    hold off
end
