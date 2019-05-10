function X = DrawMultiNorm(Miu, Sigma)

m = length(Sigma);
[P, D] = chol(Sigma);

A = P';
Z = normrnd(0, 1, m, 1);
X = Miu + A * Z;

return 