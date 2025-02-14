# XiAlphaNet
Lifespan Mapping of EEG Source Spectral Dynamics with ξ − α NET

- Authors
Ronaldo Garcia Reyes, Ariosky Areaces Gonzales, Ying Wang, Yu Jin, Ludovico Minatti and Pedro A. Valdes–Sosa

- Abstract
We formulate a new class of parametric, multivariate, and structurally informed spectral components
model of the EEG, the 𝜉 − 𝛼NET, that allows us to map the Lifespan of EEG source spectral
dynamics across a large data set of EEG cross spectrum at high spatial resolution. This approach
accurately estimates source spectral components and effective connectivity through the use of
biophysical modeling while maintaining computational efficiency, as confirmed by simulation
benchmarks. We are able to analyze source dynamics with a resolution of 8,003 voxels from the
HarMNqEEG dataset, which includes scalp EEG cross-spectrum tensors collected from 1965
subjects across 9 countries, using various devices and accounting for different age groups. Our
findings indicate that the Bayesian Model Inversion of the 𝜉 − 𝛼NET allows to map lifespan of
conduction delays that follows a U-shaped trajectory, which contrasts with independently recorded
myelin concentration measurements. Moreover, we assess the spatiotemporal distribution of spectral
components, revealing that the aperiodic or fractal component has an isotropic spatial distribution on
the cortical surface. While the generator’s spectral peak in the 𝛼 band, i.e., 𝛼-rythms, is localized on
the visual areas of the brain. Using a Zero Inflated Gaussian model, our findings indicate that the
mode frequency that characterizes the 𝛼-rythms or Peak Alpha Frequency shows an inverted
U-shaped trajectory for both hemispheres across the lifespan, and a spatial gradient of zero inflation
in PAF across the cortex that flatten the trajetory from posterior to frontal areas. We provide both the
code of the 𝜉 − 𝛼NET and the source solution of the spectral dynamics for the HarMNqEEG.