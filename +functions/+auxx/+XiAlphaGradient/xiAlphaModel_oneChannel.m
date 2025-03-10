function N_model = xiAlphaModel_oneChannel(params, freq)
% xiAlphaModel_oneChannel: Xi + Alpha model for one channel
%
%    params = [ e1, e2, e3, a1, a2, a3, a4 ]
%    freq   = frequency vector (same length as the data to fit)
%
% Returns the modeled spectrum N_model at each frequency in freq.

    % Extract parameters
    e1 = params(1);
    e2 = params(2);
    e3 = params(3);
    a1 = params(4);
    a2 = params(5);
    a3 = params(6);
    a4 = params(7);

    % Xi component: t4(omega | e1, e2, e3, 0)
    Xi = e1 ./ ( (1 + e2*(freq - 0).^2 ).^e3 );

    % Alpha component: t4(omega | a1, a2, a3, a4)
    Alpha = a1 ./ ( (1 + a2*(freq - a4).^2 ).^a3 );

    % Sum
    N_model = Xi + Alpha;
end
