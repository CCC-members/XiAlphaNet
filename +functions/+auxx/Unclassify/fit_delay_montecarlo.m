% Preprocessing
preprocessing_velocity;
load('Data/Average_Velocity_ROISpace/GPfit_Delay_Mean.mat');
delay_init = m2v(sym_matrix(corrected_delay)); % Mean delay
std = 4;

% Define the interval I for the velocity vector (v)
lb = max(delay_init - 3 * std * ones(size(delay_init)), 0); 
ub = delay_init + 3 * std * ones(size(delay_init)); % Upper bound

% Parameters
h = 100;
batch = 10;
numSamples = 10; % Number of random samples to evaluate
numVars = length(delay_init);

% Ensure lb and ub are row vectors
lb = lb(:)'; % Convert to row vector if not already
ub = ub(:)'; % Convert to row vector if not already



%%
% Preprocessing
preprocessing_velocity;
load('Data/Average_Velocity_ROISpace/GPfit_Delay_Mean.mat');
delay_init = m2v(sym_matrix(corrected_delay)); % Mean delay
std = 4;

% Define the interval I for the velocity vector (v)
lb = max(delay_init - 3 * std * ones(size(delay_init)), 0); 
ub = delay_init + 3 * std * ones(size(delay_init)); % Upper bound

% Parameters
h = 100;
batch = 20;
numSamples = 50; % Number of random samples to evaluate
numVars = length(delay_init);

% Ensure lb and ub are row vectors
lb = lb(:)'; % Convert to row vector if not already
ub = ub(:)'; % Convert to row vector if not already

% Prepare to store results
bestVs = cell(1, 10);
bestObjectiveValues = zeros(1, 10);
age_stars = linspace(10, 100, 200);

parfor a = 1:length(age_stars)
    age_star = age_stars(a);
    
    % Monte Carlo Optimization
    
    % Vectorized sampling: generate all random samples at once
    sampling_mean = (ub + lb) / 2; % Centered within the bounds
    sampling_variance = ((ub - lb) / 12).^2; % Variance proportional to the range, adjust as needed

    % Generate samples from a normal distribution
    dSamples = sampling_mean + sqrt(sampling_variance) .* randn(numSamples, numVars);

    % Ensure the samples stay within the bounds [lb, ub]
    dSamples = max(min(dSamples, ub), lb);
    

    % Initialize the best objective value and corresponding sample
    bestObjectiveValue = inf;
    bestV = [];
    
    
    % Evaluate the objective function for each sample
    for i = 1:numSamples
        objectiveValue = delay_objective_function(dSamples(i, :), age_star, h, batch, parameters);
        
        % Check if this is the best result so far
        if objectiveValue < bestObjectiveValue
            bestObjectiveValue = objectiveValue;
            bestV = dSamples(i, :);
        end
    end
   
    
    % Store the best results for this age_star
    bestVs{a} = bestV;
    bestObjectiveValues(a) = bestObjectiveValue;
end














































































% % Prepare to store results
% bestVs = cell(10, 200); % Adjusted the second dimension to match length(age_stars)
% bestObjectiveValues = zeros(10, 200);
% age_stars = linspace(10, 100, 200);
% tic
% parfor j = 1:10
%     local_bestVs = cell(1, 200); % Preallocate for local storage within parfor
%     local_bestObjectiveValues = zeros(1, 200);
%     
%     for a = 1:length(age_stars)
%         age_star = age_stars(a);
%         
%         % Monte Carlo Optimization
%         % Vectorized sampling: generate all random samples at once
%         sampling_mean = (ub + lb) / 2; % Centered within the bounds
%         sampling_variance = ((ub - lb) / 12).^2; % Variance proportional to the range, adjust as needed
%     
%         % Generate samples from a normal distribution
%         dSamples = sampling_mean + sqrt(sampling_variance) .* randn(numSamples, numVars);
%     
%         % Ensure the samples stay within the bounds [lb, ub]
%         dSamples = max(min(dSamples, ub), lb);
%     
%         % Initialize the best objective value and corresponding sample
%         bestObjectiveValue = inf;
%         bestV = [];
%         
%         % Evaluate the objective function for each sample
%         for i = 1:numSamples
%             objectiveValue = delay_objective_function(dSamples(i, :), age_star, h, batch, parameters);
%             
%             % Check if this is the best result so far
%             if objectiveValue < bestObjectiveValue
%                 bestObjectiveValue = objectiveValue;
%                 bestV = dSamples(i, :);
%             end
%         end
%         
%         % Store the best results for this age_star locally
%         local_bestVs{a} = bestV;
%         local_bestObjectiveValues(a) = bestObjectiveValue;
%     end
%     
%     % Store the results from this iteration of the parfor loop
%     bestVs(j, :) = local_bestVs;
%     bestObjectiveValues(j, :) = local_bestObjectiveValues;
% end
% toc
% 
% 

