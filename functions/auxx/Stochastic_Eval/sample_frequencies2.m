function [parameters] = sample_frequencies2(parameters)
    % This function samples frequencies randomly from the given input parameters.
    % It returns an array of sampled frequencies and their positions in the 
    % original frequency array.
    %
    % Inputs:
    %   parameters: Struct with nested structs for data and stochastic sampling
    %     Data.freq - Array of frequencies to sample from
    %     Stochastic.Nsfreq - Number of frequencies to sample
    %     Stochastic.stoch - Flag to enable stochastic sampling
    %
    % Outputs:
    %   parameters: Updated structure with sampled frequencies and their indices
    %
    % Example usage:
    %   params.Data.freq = 0.1:0.1:40;
    %   params.Stochastic.Nsfreq = 5;
    %   params.Stochastic.stoch = 1;
    %   params = sample_frequencies(params);

    % Extracting parameters
    freq = parameters.Data.freq;
    k = parameters.Stochastic.Nsfreq;
    stoch_flag = parameters.Stochastic.stoch;

    % Conditional sampling based on the stochastic flag
    if stoch_flag == 1
        % Ensure there are enough frequencies to sample from
        if length(freq) < k
            error(['Insufficient frequencies available: Requested ', num2str(k), ', but only ', num2str(length(freq)), ' available.']);
        end
        
        % Randomly sample k frequencies
        idx = randperm(length(freq), k);
        sampled_freqs = freq(idx);
        positions = idx;
    else
        % Use all frequencies if not stochastic sampling
        sampled_freqs = freq;
        positions = 1:length(freq);
    end

    % Store sampled frequencies and their indices in parameters struct
    parameters.Stochastic.Sampled.sw = [sampled_freqs(:)'; ones(1, length(sampled_freqs))];
    parameters.Stochastic.Sampled.sp = positions;
end
