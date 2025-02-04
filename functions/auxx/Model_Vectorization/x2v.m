function [e, a, sigma2] = x2v(x)
    Nr = floor((length(x)-1)/7);
    num_e = Nr * 3;  % Number of elements in e
    num_a = Nr * 4;  % Number of elements in a
    e = reshape(x(1:num_e), Nr, 3);
    a = reshape(x(num_e + 1:num_e + num_a), Nr, 4);
    sigma2 = x(num_e + num_a + 1:end);
end
