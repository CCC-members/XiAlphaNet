%% Auxx functions

addpath('auxx')
addpath('auxx/Generate_Source_Sample/')
addpath('auxx/Model_Vectorization/')
addpath('auxx/Numerical_Check/')
addpath('auxx/Optimized_Operations/')
addpath('auxx/Proximal_Operators/')
addpath('auxx/Regularization/')
addpath('auxx/Stochastic_Eval/')
addpath('auxx/T_Operator/')
addpath('auxx/XiAlpha_Gradient/')
addpath('auxx/Data_Preprosessing/')
addpath('auxx/Gaussian Process Regression/')

%% Data 

addpath('Data')
addpath('Data/Atlas_Anatomical/')
addpath('Data/Average_Conn&Tract_VoxelSpace/')
addpath('Data/Avarege_Conn&Tract_ROISpace/')
addpath('Data/Lead_Field/')
addpath('Data/Model_Parameters/')
addpath('Data/Scalp_Density_Matrix/')
addpath('Data/Average_Velocity_ROISpace/')

%% Fig 
%addpath('fig/')
%% Objective function, Gradient & Proximal Operator 

addpath('Function_Grandient_Prox')
addpath('Function_Grandient_Prox/Eval_Gradient/')
addpath('Function_Grandient_Prox/Eval_Prox_Operator/')
addpath('Function_Grandient_Prox/Eval_Regularization_Term/')
addpath('Function_Grandient_Prox/Eval_Smooth_Term/')

%% Source Crosspectra

addpath('Generate_Source_Conn')

%% Stocastich FISTA

addpath('Stochastic_FISTA')

%% Visualization 

addpath('Visualization/')
addpath('Visualization/BC-VARETA_Toolbox_viewer/')
addpath('Visualization/ESI_plot/')
addpath('Visualization/DataViz-master/daviolinplot/')
