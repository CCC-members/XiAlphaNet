function [p] = prox(x, c)
import functions.*
import functions.auxx.*
import functions.auxx.ModelVectorization.*
import functions.auxx.ProximalOperators.*

% Split x into parts e, a, s2 using x2v function
[e, a, s2] = x2v(real(x));

% Retrieve the threshold parameter

% Apply ridge and group lasso proximal functions
% Ensure s2 is above a minimum value
s2 = abs(max(prox_ridge(s2, c(1)), 10^(-3)));
% Ensure e and a are non-negative after applying proximal functions
e(:,1:3) = abs(max(prox_glasso(e(:,1:3), c(2)), 0));
a(:,1:3) = abs(max(prox_glasso(a(:,1:3), c(3)), 0));

% Identify non-zero entries in a and e
nonzero_positions_a = a(:, 1) ~= 0;
nonzero_positions_e = e(:, 1) ~= 0;

%     % Apply conditional modifications based on non-zero positions in a
if any(nonzero_positions_a)
    a(nonzero_positions_a, 4) = max(min(a(nonzero_positions_a, 4), 13), 7);
    a(nonzero_positions_a, 2) = max(a(nonzero_positions_a, 2), 0);
    a(nonzero_positions_a, 3) = max(a(nonzero_positions_a, 3),  0);
    % a(nonzero_positions_a, 1) = min(a(nonzero_positions_a,1),e(nonzero_positions_a,1)/5);
end
%
%     % Apply conditional modifications based on non-zero positions in e
if any(nonzero_positions_e)
    e(nonzero_positions_e, 2) = max(e(nonzero_positions_e, 2),  0);
    e(nonzero_positions_e, 3) = max(e(nonzero_positions_e, 3), 0);
end

% Identify zero entries in a and e
zero_positions_a = a(:, 1) == 0;
zero_positions_e = e(:, 1) == 0;

% Set columns 2 to 4 to zero for zero positions in a
if any(zero_positions_a)
    a(zero_positions_a, 2:4) = 0;
end

% Set columns 2 and 3 to zero for zero positions in e
if any(zero_positions_e)
    e(zero_positions_e, 2:3) = 0;
end

% Identify zero entries in a and e
zero_positions_a2 = a(:, 4) == 0;

% Set columns 2 to 4 to zero for zero positions in a
if any(zero_positions_a2)
    a(zero_positions_a2, 1) = 0;
end


% Reassemble the parts e, a, s2 back into p using v2x function
p = v2x(e, a, s2);
end
