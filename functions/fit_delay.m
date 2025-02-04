%% Preprocessing
preprocessing_velocity;
load('Data/Average_Velocity_ROISpace/GPfit_Delay_Mean.mat');
delay_init = m2v(sym_matrix(corrected_delay)); % Mean delay
std = 4;

%% Define the interval I for the velocity vector (v)
lb = max(delay_init - 3 * std * ones(size(delay_init)), 0); 
ub = delay_init + 3 * std * ones(size(delay_init)); % No need for max, since it should be non-negative by construction

%% Create the optimization variables
% Initialize an empty array to store the optimization variables
optVars = [];


% Preallocate the cell array for optimizable variables
numVars = length(delay_init);
optVars = cell(1, numVars); % Preallocate a cell array

% Fill in the cell array
for i = 1:numVars
    i
    % Define an optimizable variable for each element in the vector
    varName = ['d', num2str(i)];
    optVars{i} = optimizableVariable(varName, [lb(i), ub(i)], 'Type', 'real');
end
1
% Convert the cell array to a standard array for use in bayesopt
optVars = [optVars{:}];
2
%% Set up the objective function
age_star = 10;
h = 100;
batch = 1;
objectiveFunction = @(d) delay_objective_function(struct2array(d), age_star, h, batch, parameters);

%% Set up the Bayesian Optimization options
results = bayesopt(objectiveFunction, optVars, ...
    'MaxObjectiveEvaluations', 30, ... % Adjust as necessary
    'IsObjectiveDeterministic', true, ...
    'AcquisitionFunctionName', 'expected-improvement-plus', ...
    'Verbose', 1); % Set verbosity level to 1 (basic information)

%% Extract the best value
bestV = table2array(results.XAtMinObjective);
