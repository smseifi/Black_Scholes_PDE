%% 1. Finite Difference, Constant Timesteps, European 

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

r = 2e-2;
T = 75e-2;
K = 5e+1;
S_zero = 5e+1;

localSigma = @(S) 85e-2./sqrt(S);

n = 25; delTau = T/n;

load coarse_mesh.mat
V_N = chooserPayoff(S_NVect, K);

sigma = localSigma(S_NVect);

L = discreteBls(S_NVect, r, sigma);
Value = fdSolve(V_N, L, delTau, n);

nn = 5;
ATM = zeros(nn+1, 1);
ATM(1, 1) = Value(S_NVect==S_zero, 1);
for i = 1:nn
    n =2*n; delTau = delTau/2;
    m = size(S_NVect, 1);
    S = zeros(2*m-1, 1);
    S(1:2:end, 1) = S_NVect;
    S(2:2:end, 1) = (S_NVect(2:end, 1) + S_NVect(1:end-1, 1))/2;
    V_N = chooserPayoff(S, K);
    sigma = localSigma(S);
    L = discreteBls(S, r, sigma);
    Value = fdSolve(V_N, L,delTau, n);
    ATM(i+1, 1) = Value(S==S_zero, 1);
    S_NVect = S;
end