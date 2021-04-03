function soln = implicit(V, L, dt, n)
% Description of implicit goes here
%   Detailed description goes here

m = size(V, 1);
soln = zeros(m, n+1);  soln(:, n+1) = V;

for i = n:-1:1
    soln(:, i) = (eye(m, m) - dt*L) \ soln(:, i+1);
end
end

