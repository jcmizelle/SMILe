# SMILe EEG Toolbox (MATLAB)

MATLAB code from the Sensorimotor Integration Lab (SMILe) for EEG preprocessing, spectral and time-frequency analysis, functional connectivity, graph theory, surrogate data, and permutation statistics. Written by J. Christopher Mizelle, PhD, for lab studies of brain injury, aging, and sensorimotor control.

About 17,000 lines across 70 functions and scripts. Most functions follow one pattern: data in as `chan x time x trial` arrays (or an EEGLAB `EEG` struct), options in a `cfg` struct with documented defaults, and an input/output header.

## Modules

| Folder | What it does |
|---|---|
| `SMILe_Preprocessing` | Config-driven pipeline on top of EEGLAB: import (CNT/CDT/SET), resampling, channel handling, high/low-pass filtering, line-noise removal, ASR artifact cleaning, interpolation, re-referencing (average/REST), epoching, ICA (AMICA, SOBI, JADE, Picard) with ICLabel/MARA classification, baseline correction, dipole fitting, surface Laplacian (CSD), and time-frequency output |
| `SMILe_Connectivity_and_Graph` | Cross-spectra from Fourier, wavelet, and Hilbert/Gabor decompositions; connectivity metrics (PLV, iPLV, ciPLV, partial PLV, instantaneous PLV, coherence, imaginary/lagged coherence, PSI, frontal alpha asymmetry); permutation z-thresholding of adjacency matrices; surrogate testing; minimum spanning trees (Prim, Kruskal, Euclidean); graph-measure wrappers; network visualizations (circle, line, arrow, dot, property, and hub topoplots) |
| `SMILe_Frequency_Utilities` | Convolution, FIR, and wavelet band filters; amplitude envelopes; simulated connectivity trials |
| `SMILe_Surrogates` | Surrogate data generation: phase randomization (single, epoch-wise), AAFT, iAAFT, time permutation, shuffling, Box-Cox |
| `SMILe_Statistics`, `SMILe_Resampling` | Permutation z-thresholding for time-voltage and time-frequency data, bootstrap confidence intervals, surrogate simulations |
| `SMILe_Factorial_Mass_Univar_Test` | Driver for factorial mass-univariate ERP analysis (FMUT/MUT) |
| `SMILe_Plotting` | Topographic maps, hub topoplots, and video rendering of time-resolved maps |
| `SMILe_Power_Spectrum`, `SMILe_SourceImaging`, `SMILe_Infinity_Reference`, `SMILe_Processing_Utilities`, `SMILe_Workflows` | Power spectra, distributed source modeling, REST reference, baseline removal, and an end-to-end connectivity + graph workflow |

## Dependencies (not included)
- [EEGLAB](https://sccn.ucsd.edu/eeglab/) with clean_rawdata, cleanline, ICLabel, AMICA, picard, MARA, SASICA, dipfit, REST
- [FieldTrip](https://www.fieldtriptoolbox.org/) (layouts and some plotting)
- [Brain Connectivity Toolbox](https://sites.google.com/site/bctnet/) (graph measures)
- [CSD Toolbox](https://psychophysiology.cpmc.columbia.edu/software/CSDtoolbox/) (surface Laplacian)
- [FMUT](https://github.com/ericcfields/FMUT) / [Mass Univariate ERP Toolbox](https://github.com/dmgroppe/Mass_Univariate_ERP_Toolbox)

Third-party code that lived alongside these functions in the lab copy has been left out; use the original sources above.

## Notes
- File paths are placeholders (`path/to/...`). Point them at your data and electrode templates.
- Some scripts are lab working scripts rather than finished functions and expect EEGLAB's base workspace.

## License
GPL-3.0. See LICENSE.
