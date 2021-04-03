%% 1. Finite Difference, Constant Timesteps, European 
% ----------------------------------------------------------------------- %
c = input(['Choose one of the following time-stepping schemes' ... 
'\n 1: Implicit time-stepping \n 2: Crank--Nicolson time-stepping' ...    
'\n 3: Crank--Nicolson--Rannacher time-stepping \n']);

switch c
    case 1
        fdSolve = @(V, L, dt, n) implicit(V, L, dt, n);
    case 2
        fdSolve = @(V, L, dt, n) cn(V, L, dt, n);
    case 3
        fdSolve = @(V, L, dt, n) cnr(V, L, dt, n);
end
% ----------------------------------------------------------------------- %

r = 2e-2;
T = 75e-2;
K = 5e+1;
S_zero = 5e+1;
localSigma = @(S) 85e-2./sqrt(S);
% ----------------------------------------------------------------------- %

load coarse_mesh.mat

nn = 5;
nVect = 25 * 2.^(0:nn).';
deltVect = T./nVect; 

mVect = zeros(nn+1, 1); mVect(1, 1) = size(S_NVect, 1);
ATM = zeros(nn+1, 1);

dt = deltVect(1, 1);    n = nVect(1, 1);
V_N = chooserPayoff(S_NVect, K);
sigma = localSigma(S_NVect);
L = discreteBls(S_NVect, r, sigma);
Value = fdSolve(V_N, L, dt, n);

ATM(1, 1) = Value(S_NVect==S_zero, 1);
for i = 2:nn+1
    n =nVect(i, 1); dt = deltVect(i, 1); 
    m = size(S_NVect, 1);
    S = zeros(2*m-1, 1);
    mVect(i, 1) = 2*m-1;
    S(1:2:end, 1) = S_NVect;
    S(2:2:end, 1) = (S_NVect(2:end, 1) + S_NVect(1:end-1, 1))/2;
    V_N = chooserPayoff(S, K);
    sigma = localSigma(S);
    L = discreteBls(S, r, sigma);
    Value = fdSolve(V_N, L,dt, n);
    ATM(i, 1) = Value(S==S_zero, 1);
    S_NVect = S;
end
% ----------------------------------------------------------------------- %

change = diff(ATM); ratio = change(1:nn-1) ./ change(2:nn);
change = [nan; change]; ratio = [nan*ones(2, 1); ratio];

tb1 = table(mVect, nVect, ATM, change, ratio, 'VariableNames', ...
    {'Nodes', 'Time Stamps', 'ATM Chooser Value', 'Change', 'Ratio'});
disp(tb1)