function soln = cnr(V, L, dt, n)
% Description of cnr goes here
%   Detailed description goes here

m = size(V, 1);
soln = zeros(m, n+1);  soln(:, n+1) = V;

for i = 1:min(2, n)
    soln(:, n+1-i) = (eye(m, m) - dt*L) \ soln(:, n+2-i);
end
for i = min(3, n):n
    soln(:, n+1-i) = (eye(m, m) - 5e-1*dt*L) \ ((eye(m, m) + ...
        5e-1*dt*L) * soln(:, n+2-i));
end
end
