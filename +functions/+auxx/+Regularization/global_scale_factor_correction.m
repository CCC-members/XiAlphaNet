function [S,scale] = global_scale_factor_correction(S);
scale = log(max(abs(S(:))));
factor = 20-scale;
S = S*(exp(factor));
scale = exp(factor);
end