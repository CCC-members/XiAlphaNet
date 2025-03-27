
# Xi-AlphaNET

**Version**: 1.0  
**Release Date**: February 2025  
**Authors**: Ronaldo Garcia Reyes, Ariosky Areaces Gonzales, Pedro A. Valdes–Sosa, Xi-AlphaNET Research Team

---

## Overview

Xi-AlphaNET is a parametric, multivariate model developed to map the lifespan of EEG source spectral dynamics across large datasets with high spatial resolution. This tool uses biophysical modeling to estimate source spectral components and effective connectivity while maintaining computational efficiency. The model integrates aperiodic and periodic spectral components, specifically the alpha rhythm (α), and characterizes their spatiotemporal distributions across cortical regions. By analyzing EEG data from a global cohort, the tool allows researchers to explore age-related changes in brain activity, including the unique U-shaped trajectory of conduction delays and the isotropic spatial distribution of aperiodic components.

---

## Output Files Overview

After processing each subject, Xi-AlphaNET generates the following files in the output directory:

| File Name                     | Description                                                                                                                                           |
|-------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------|
| `Alpha_estimate.mat`           | Contains voxel-wise alpha parameters: Power, Width, Exponent, and Peak Alpha Frequency.                                                                |
| `Xi_estimate.mat`              | Contains voxel-wise aperiodic "Xi" parameters, structured similarly to Alpha.                                                                         |
| `Conn_Matrix.mat`              | Latent anatomical connectivity matrix (Region x Region), based on the HCP-MMP1 atlas unless customized.                                                |
| `Delay_Matrix.mat`             | Estimated neural conduction delays between regions (Region x Region).                                                                                 |
| `Mod_Weights.mat`              | Metadata vector (2 components) regarding model regularization and fitting quality.                                                                   |
| `Source_Activations_Full.mat`  | Total source activations (Voxel x Frequency bins). Equal to Alpha + Xi activations.                                                                  |
| `Source_Activations_Alpha.mat` | Source activation attributable to alpha peaks.                                                                                                        |
| `Source_Activations_Xi.mat`   | Source activation attributable to aperiodic (Xi) background activity.                                                                                 |
| `Source_PSD.mat`               | Estimated Power Spectral Density (ROI x Frequency). ROI-level resolution unless full resolution is selected.                                           |

---

## Step-by-Step User Workflow

### 1. Data Preparation
- Supported formats: Raw EEG (.PLG, .EDF, .MAT) or preprocessed cross-spectrum .mat files.
- Organize your data into folders: one subfolder per subject.
- Run `Update_data.m` from the `/utils/` folder to structure the dataset:
  - Define `input_folder` (your current data location).
  - Define `output_folder` (where structured folders will be saved).

### 2. Generate Participant Info
- Run `get_participants_info.m`.
- Outputs a `Participant.json` file with demographic information (Name, Age, Sex, etc.).
- Input path should point to the folder structured in Step 1.

### 3. Launch the Application
- Launch via:
  - MATLAB console: `xialphanet`
  - Windows executable: double-click the `.exe`

### 4. App Interface Guide

#### 4.1 Load Saved Data (Optional)
- Upon launch, Xi-AlphaNET checks for previously saved analyses.
- Select:
  - `Load`: Load existing results into the interface.
  - `Close`: Start fresh with a new dataset.

### 5. Configure Settings

#### General Tab
- **Dataset Name**: Set project name.
- **Description**: Optional notes about the dataset.
- **Input Path**: Folder containing EEG or cross-spectrum files (from Step 1).
- **Reference File**: Select an example `.PLG`, `.EDF`, or `.MAT` to define file structure.
- **Participant File**: Load `Participant.json` created in Step 2.
- **Data Type**: Choose between EEG Signal or Cross-Spectrum.

#### Output & Computation Paths
- **Output Path**: Where results will be stored.
- **Tmp Path**: `local` by default. But it can be changed to the location of the Folder that contains the TensorField if the user decides to download these values and use the default anatomical input. 
- **Parallel Computation**: Enabled by default. Can be toggled.
- **Save Mode**: Choose to store extended outputs such as full cross-spectra.

### 6. Anatomy Tab
- Uses the HCP-MMP1 template by default.
- Option to disable the default and provide custom anatomical inputs:
  - Cortex surface
  - LeadField
  - Channel file
  - AnatomicalConnectivity
  - Neuro Tract Length
  - Conduction Delay Matrix
- **Important**: All anatomical data must be in the same spatial reference system and match in dimension.

### 7. Model Parameters Tab

| Field                      | Description                                                                                                                                 |
|----------------------------|---------------------------------------------------------------------------------------------------------------------------------------------|
| **Use Default Tensor Field**| Downloads ~30GB precomputed transfer functions (recommended for large datasets with default anatomy).                                         |
| **Number of Frequencies**   | Choose the number of frequency bins. E.g., 47 bins ≈ 20 Hz for HarmNqEEG.                                                                   |
| **Bayes Iteration Delays**  | Controls optimization depth for delay/connectivity search. Recommended: ≥30 (ideal: 100–200).                                               |
| **Nrand 1 & Nrand 2**       | Number of random initializations: for regularization search and final fitting, respectively.                                               |
| **Bayes Iteration 1 & 2**   | Control search for sparsity regularization before and after model fitting.                                                                  |
| **Stochastic Evaluation**   | Toggles batch-wise stochastic FISTA optimization. Recommended ON.                                                                           |
| **Age Groups**              | Define or edit age-based participant groups (optional).                                                                                     |

### 8. Running the Model
- Click **Apply** to save settings.
- Then click **Run** (top-left).

**Note**: Runtime can vary:
- One subject ≈ 0.6 minutes on a high-performance machine (256 GB RAM, 52 cores) using default Tensor Field.
- Longer if the default Tensor Field is not used around 1.5 - 2 minutes.

### 9. Viewing Results
- After processing is complete, the dataset icon will appear under the `Datasets` field.
- Select a dataset to expand:
  - **Left**: List of all subjects.
  - **Right**: List of visualizations (e.g., Alpha Power Map).

### 10. Visualization Options
- Click on a subject to view:
  - Alpha/Xi power distribution on the cortex
  - Delay and connectivity matrices
  - Source PSD across frequencies
  - Spatial maps of alpha peak frequency, width, exponent

---

## Tips & Troubleshooting

| Issue                     | Solution                                                                                  |
|---------------------------|-------------------------------------------------------------------------------------------|
| **No subjects detected**   | Check input path and reference file. Must point to correct folder and format.             |
| **Missing outputs**        | Verify you’ve selected the proper data type (EEG vs Cross) and saved correctly.            |
| **High memory usage**      | Use the precomputed Tensor Field and enable parallel computation.                         |
| **Slow analysis**          | Reduce Bayes Iterations,Nrand1 and Nrand2 for faster but less precise results. Additionally if you use the defaults anatomical templates, consider download the  Tensor Field for speed up calculations.                              |

---

## Contact & Resources
- **Documentation**: https://www.researchgate.net/publication/389324537_Lifespan_Mapping_of_EEG_Source_Spectral_Dynamics_with_x_-_aNET
- **Tensor Field Download**: Provided on the website (~30 GB)
- **Support Contact**: ronaldo@neuroinformatics-collaboratory.org
