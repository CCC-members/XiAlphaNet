
# Xi-AlphaNET

**Version**: 5.0  
**Release Date**: November 2025  
**Authors**: Ronald Garcia Reyes, Ariosky Areaces Gonzales, Pedro A. Valdes–Sosa, Xi-AlphaNET Research Team

---

## Overview

Xi-AlphaNET is a parametric, multivariate generative framework designed to map lifespan trajectories of EEG source spectral dynamics at high spatial resolution across large populations. Grounded in biophysical modeling, the framework jointly estimates source-level spectral parameters, effective connectivity, and interareal conduction delays while remaining computationally scalable. In its current formulation, Xi-AlphaNET explicitly models aperiodic (ξ) activity and the alpha rhythm (α) as anatomically constrained stochastic processes, enabling characterization of their distinct spatiotemporal organization across the cortex. Applied to large, harmonized EEG cohorts, the model reveals systematic age-dependent effects, including nonlinear inverted-U trajectories of spectral parameters and global conduction delays, as well as an isotropic spatial organization of aperiodic activity. Future releases will extend this framework to incorporate additional spectral components, providing a more comprehensive generative description of cortical electrophysiological dynamics across the lifespan.


## Output Files Overview

After processing each subject, Xi-AlphaNET produces a standardized set of output files stored in the designated results directory. These files contain source-level spectral estimates, connectivity measures, and model-specific parameters:

| File Name                      | Description                                                                                                                                                                                                                                     |
| ------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `Alpha_estimate.mat`           | Voxel-wise estimates of the alpha-rhythm spectral parameters, including amplitude Aα, bandwidth Bα, spectral exponent Eα, and peak alpha frequency Fα. Each parameter is defined at the cortical vertex level and has dimension N_vertices × 1. |
| `Xi_estimate.mat`              | Voxel-wise estimates of the aperiodic (ξ) spectral parameters, including amplitude Aξ, bandwidth Bξ, and spectral exponent Eξ. Each parameter is defined at the cortical vertex level and has dimension N_vertices × 1.                         |
| `Conn_Matrix.mat`              | Estimated latent anatomical connectivity matrix (Region × Region), defined on the HCP-MMP1 parcellation unless a custom atlas is provided.                                                                                                      |
| `Delay_Matrix.mat`             | Estimated inter-regional neural conduction delays (Region × Region).                                                                                                                                                                            |
| `Mod_Weights.mat`              | Estimated regularization parameters of the Xi-AlphaNET generative model, specifically the weights wC and wD controlling the contribution of anatomical connectivity and conduction delays, respectively.                                        |
| `Source_Activations_Full.mat`  | Total source-level activations across frequencies (N_vertices × Frequency), corresponding to the combined contribution of the α and ξ processes.                                                                                                |
| `Source_Activations_Alpha.mat` | Source-level activations attributable specifically to the alpha rhythmic process (N_vertices × Frequency).                                                                                                                                      |
| `Source_Activations_Xi.mat`    | Source-level activations attributable to the aperiodic (ξ) background process (N_vertices × Frequency).                                                                                                                                         |
| `Source_PSD.mat`               | Estimated power spectral density at ROI or vertex resolution, depending on the selected configuration.                                                                                                                                          |


## Step-by-Step User Workflow

### 1. Data Preparation
- Supported input formats include raw EEG recordings (`.PLG`, `.EDF`, `.MAT`, `.set`) and, optionally, precomputed cross-spectral data stored in `.mat` format. Plain-text (`.txt`) input will be supported in future releases.
- For optimal results, EEG preprocessing should be performed **on raw EEG signals** prior to any spectral estimation. Artifact rejection, bad-channel identification, re-referencing, filtering, and epoch selection should be carried out using established EEGLAB functions.
- The use of precomputed cross-spectra is supported for advanced users; however, preprocessing at the cross-spectral level is not recommended, as it limits artifact control and channel-level quality assessment.
- Organize the dataset such that each subject is contained in a separate subfolder.
- Run `Update_data.m` from the `/utils/` directory to standardize the dataset structure:
  - Define `input_folder` as the directory containing the original EEG data.
  - Define `output_folder` as the directory where the structured dataset will be created.

### 2. Generate Participant Information
- Run `get_participants_info.m`.
- This step generates a `Participant.json` file containing participant metadata (e.g., age, sex).
- The input path must correspond to the structured dataset created in Step 1.

### 3. Launch the Application
- Launch Xi-AlphaNET using one of the following methods:
  - MATLAB command window: `xialphanet`

### 4. Application Interface

#### 4.1 Load Saved Data (Optional)
- Upon startup, Xi-AlphaNET automatically checks for previously saved analyses.
- Select:
  - **Load** to open existing results in the interface.
  - **Close** to start a new analysis with a fresh dataset.

### 5. Configure Settings

#### General Tab
- **Dataset Name**: Specify a project identifier.
- **Description**: Optional textual description of the dataset.
- **Input Path**: Directory containing the EEG or cross-spectral data (from Step 1).
- **Reference File**: Select an example EEG file (`.PLG`, `.EDF`, `.MAT`, or `.set`) to define the expected data structure.
- **Participant File**: Load the `Participant.json` file generated in Step 2.
- **Data Type**: Select either raw EEG signals or cross-spectral data.

#### Output & Computation Paths
- **Output Path**: Directory where all results will be stored.
- **Temporary Path**: Location used for intermediate files during computation (local by default).
- **Parallel Computation**: Enabled by default and recommended for large datasets; can be disabled if required.
- **Save Mode**: Option to store extended outputs, such as intermediate spectral quantities.

### 6. Anatomy Tab
- By default, Xi-AlphaNET uses the **HCP-MMP1 cortical parcellation**, which defines the spatial resolution of the anatomical connectivity and conduction delay matrices.
- Users may optionally disable the default anatomy and provide custom anatomical inputs:
  - Cortical surface
  - Lead field
  - Channel definition file
  - Anatomical connectivity matrix
  - Tract length information
  - Conduction delay matrix
- **Important**: All anatomical inputs must be defined in the same spatial reference frame and be dimensionally consistent.

### 7. Model Parameters Tab

| Field | Description |
|------|-------------|
| **Number of Frequencies** | Number of frequency bins used for spectral estimation (e.g., 47 bins ≈ 20 Hz for HarMNqEEG data). |
| **Bayes Iteration Delays** | Controls the optimization depth for estimating effective connectivity and conduction delays. Recommended ≥ 30 (typical range: 100–200). |
| **Nrand 1 & Nrand 2** | Number of random initializations used during regularization search and final model fitting, respectively. |
| **Bayes Iteration 1 & 2** | Controls sparsity regularization before and after model fitting. |
| **Stochastic Evaluation** | Enables batch-wise stochastic optimization; recommended for large datasets. |
| **Age Groups** | Define or edit age-based participant groupings (optional). |

### 8. Running the Model
- Click **Apply** to save the current configuration.
- Click **Run** to start model estimation.

**Note**: Runtime depends on dataset size, parameter configuration, and available computational resources. The software scales efficiently on multi-core systems using parallel computation.

### 9. Viewing Results
- Upon completion, the processed dataset appears in the **Datasets** panel.
- Expanding a dataset displays:
  - **Left panel**: List of processed subjects.
  - **Right panel**: Available visualizations (e.g., alpha power maps).

### 10. Visualization Options
- Selecting a subject enables visualization of:
  - Cortical distributions of α and ξ power
  - Effective connectivity and conduction delay matrices
  - Source-level power spectral densities
  - Spatial maps of spectral parameters, including peak alpha frequency, bandwidth, and spectral exponent

---


## Contact & Resources
- **Primary reference**:  
 [ https://www.researchgate.net/publication/389324537_Lifespan_Mapping_of_EEG_Source_Spectral_Dynamics_with_x_-_aNET](https://www.biorxiv.org/content/biorxiv/early/2025/08/01/2025.02.21.639413.full.pdf)
- **Support contact**:  
  ronaldo@neuroinformatics-collaboratory.org


