function L = discreteBls(S, r, sigma)
% Description of disBls goes here
%   Detailed description goes here

[alpha, beta] = prep(S, r, sigma);

m = size(S, 1);
d = zeros(m, 1);    d(1, 1) = -r;    d(2:end-1) = -(alpha + beta + r);
bd = [alpha; zeros(2, 1)];
ad = [zeros(2, 1); beta];

L = spdiags([bd, d, ad], -1:1, m, m);
end

function [alpha, beta] = prep(S, r, sigma)
% Description of prep goes here
%   Detailed description goes here

m = size(S, 1);
inner = 2:m-1;    left = 1:m-2; right = 3:m;

num = (sigma(inner) .* S(inner)).^2;
alphDenom = (S(inner) - S(left)) .* (S(right) - S(left));
betaDenom = (S(right) - S(inner)) .* (S(right) - S(left));
centExt = r * S(inner) ./ (S(right) - S(left));
forwExt = r * S(inner) ./ (S(right) - S(inner));

alpha = num ./ alphDenom - centExt;
beta = num ./ betaDenom + centExt;

alphaNegIdx = alpha < 0;
if any(alphaNegIdx)
    alphaNegIdx = alpha < 0;
    alpha(alphaNegIdx) = alpha(alphaNegIdx) + centExt(alphaNegIdx);
    beta(alphaNegIdx) = beta(alphaNegIdx) - centExt(alphaNegIdx) + ...
        forwExt(alphaNegIdx);
end
end