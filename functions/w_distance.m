function [d] = w_distance(Sr,Se)
d = log(abs(det(inv(Se)))) -  real(trace(Sr\Se));
end
