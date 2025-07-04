function C = coherence(S)
% Assuming S is the cross-spectrum tensor of size [N_channels, N_channels, N_f]
% S(i, j, :) gives the cross-spectrum between channels i and j at all frequencies

% Number of channels and frequencies
[N_channels, ~, N_f,N_dummy] = size(S);
S = reshape(S,[N_channels,N_channels,N_f]);

% Initialize the coherence tensor
C = zeros(N_channels, N_channels, N_f);

% Loop over all pairs of channels
for i = 1:N_channels
    for j = 1:N_channels
        % Compute the auto-spectra (diagonal elements)
        Sxx = abs(S(i, i, :)).^2;  % Auto-spectrum of channel i
        Syy = abs(S(j, j, :)).^2;  % Auto-spectrum of channel j
        Sxy = abs(S(i, j, :)).^2;  % Cross-spectrum between channels i and j

        % Compute the coherence for each frequency
        C(i, j, :) = Sxy ./ (Sxx .* Syy);
    end
end
