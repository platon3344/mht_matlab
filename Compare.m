
noPlot = 0;

[state, meas, T, nStep, numTarget] = StraightFour(noPlot);

estmMHT = MHT(state, meas, T, nStep, numTarget, noPlot);
fprintf('MHT DONE\n');
