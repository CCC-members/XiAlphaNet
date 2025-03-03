function x = v2x(e, a, sigma2)
    e_flat = e(:);
    a_flat = a(:);
    x = [e_flat; a_flat; sigma2];
end
