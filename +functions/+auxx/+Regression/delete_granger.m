% Paths
subj_folder   = "/Users/ronald/Downloads/xialphanet_Solutions/Post/AMABLE_GARRIDO";
params_folder = "/Users/ronald/Downloads/xialphanet_Solutions/Post/structural";


% Compute effective connectivity
[R,~]=functions.auxx.ModelVectorization.roi_average_operators(Cortex,10);
[Cxi, Calpha] = functions.auxx.CrossSpectrum.build_effective_connectivity(subj_folder, params_folder,R);
