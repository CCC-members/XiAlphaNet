function [Nsfreq,sw,sp] = sample_frequencies(freq,index_stoch,Nsfreq)
    % This function samples frequencies randomly from the given input parameters
    % using stratified sampling if enabled. It returns an array of sampled frequencies and their
    % positions in the original frequency array.
    %
    % Inputs:
    %   parameters: Struct with two main fields
    %     Data.freq - Array of frequencies to sample from
    %     Stochastic.Nsfreq - Number of frequencies to sample
    %     Stochastic.stoch - Flag to enable/disable stochastic sampling
    %
    % Outputs:
    %   The function updates the parameters structure by adding sampled frequencies
    %   and their positions.

    % freq = parameters.Data.freq;

    if index_stoch
        k =Nsfreq;

        % Ensure there are enough frequencies to sample
        if numel(freq) < k
            error(['Insufficient frequencies available. Requested: ', num2str(k), ', Available: ', num2str(numel(freq))]);
        end
        
        % Randomly select k unique indices using stratified sampling
        numFreq = numel(freq);
        edges = linspace(1, numFreq + 1, k + 1);
        idx = zeros(1, k);  % Initialize index array

        for i = 1:k
            % Define the range of each stratum
            stratumStart = floor(edges(i));
            stratumEnd = floor(edges(i + 1)) - 1;

            % Randomly sample one frequency from each stratum
            idx(i) = randi([stratumStart, stratumEnd]);
        end

        sampled_freqs = freq(idx);
        positions = idx;
        N_numerator = numel(freq);

        % Store the sampled frequencies and their indices
        sw = [sampled_freqs; N_numerator*ones(1, k)];
        sp = positions;
    else
        % If stochastic sampling is not enabled, return all frequencies with an identity marker
        sampled_freqs = freq;
        positions = 1:numel(sampled_freqs);
        N_numerator = 1;  % Marker is set to 1 when no stochastic sampling
        Nsfreq  = 1;
        sw = [sampled_freqs; N_numerator*ones(1, numel(sampled_freqs))];
        sp = positions;
    end
end
