function V = chooserPayoff(S, K)
% Description of chooserPayoff goes here
%   Detailed description goes here

V = max(S-K, K-S);

end

