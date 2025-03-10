% Set the path to the directory containing the .mat files
mat_dir = 'C:\Users\pedro\Desktop\Phy-SCE\Data\Scalp_Density_Matrix\Control'; % Replace with your directory

% Get a list of all .mat files in the directory
files = dir(fullfile(mat_dir, '*.mat'));

% Initialize variables to store the values of Gs
Gs_values = [];

% Loop through each .mat file and process it
for file_idx = 1:length(files)
    % Load the .mat file
    load(fullfile(mat_dir, files(file_idx).name));

    % Assuming 'data_struct' is loaded from the .mat file, extract the necessary data
    Nw = 47;
    S = data_struct.CrossM(:,:,1:Nw);
    
    % Initialize the variable for k
    k = 0;
    
    % Compute k based on the provided formula
    for j = 1:Nw
        k = k + real(trace(diag(log(diag(S(:,:,j))))));
    end
    
    % Compute Gs for this .mat file
    Nc = 19;
    k = k / (Nc * Nw);
    Gs = exp(k);
    
    % Store the Gs value for this file
    Gs_values = [Gs_values, Gs];
end

% Compute the mean and standard deviation of Gs across all files
Gs_mean = mean(Gs_values);
Gs_std = std(Gs_values);

% Display the results
disp(['Mean Gs: ', num2str(Gs_mean)]);
disp(['Std Gs: ', num2str(Gs_std)]);
