%--------------------------------------------------------------------------
%
% This collection of code works on preprocessed data (typically coming out
% of the SMILe_Preprocessing pipeline that uses EEGLAB for most of the
% heavy lifting. These scripts do not rely on any external toolboxes.
%
% All functions will return a description, and also input and output
% information with the 'help' call, for example: help SMILe_dswPLI
% 
% Once cleaned and epoched, the data are then ready for frequency
% decomposition, which is the purpose of the scripts enclosed in the
% various subfolders beneath this one.
%
%--------------------------------------------------------------------------
% Cross-Spectral Analysis
%--------------------------------------------------------------------------
%
% First, the data are brought in and the cross spectrum is estimated:
%
%   - SMILe_CS_Fourier 
%   - SMILe_TF_Wavelet 
% 
% There are a number of parameters to set, and each are briefly explained 
% in the body of the code.
%
%--------------------------------------------------------------------------
% Connectivity Analysis
%--------------------------------------------------------------------------
%
% Then, the cross spectral estimates generated through these functions
% can be used to calculate a number of different functional connectivity 
% measures, based on both frequency and phase:
%
%   - SMILe_COH: returns complex cohrerency and magnitude squared
%       coherence; requires the complex cross spectrum with dimensions
%       (nchan x nchan x [nfreq x nepochs])
%
%   - SMILe_ciCOH: returns imaginary coherence and the corrected form of 
%       imaginary coherence; requires the complex cross spectrum with 
%       dimensions (nchan x nchan x [nfreq x nepochs])
%
%   - SMILe_laggedCOH: returns lagged coherence; requires the complex cross 
%       spectrum with dimensions (nchan x nchan)
%
%   - SMILe_PLV: returns the phase locking value; requires the complex 
%       cross spectrum with dimensions (nchan x nchan x nepochs)
%
%   - SMILe_instantPLV: returns the instantaneous phase locking value as a 
%       time series; requires bandpass filtered time series data as input
%
%   - SMILe_partialPLV: returns the partial phase locking value; ; requires 
%       the complex cross spectrum with dimensions (nchan x nchan)
%
%   - SMILe_PLI: returns the phase lag index; requires the complex cross 
%       spectrum with dimensions (nchan x nchan x nepochs)
%
%   - SMILe_wPLI: returns the weighted phase lag index; requires the 
%       complex cross spectrum with dimensions (nchan x nchan x nepochs)
%
%   - SMILe_dsPLI: returns the debiased squared phase lag index; requires 
%       the complex cross spectrum with dimensions 
%       (nchan x nchan x nepochs)
%
%   - SMILe_dswPLI: returns the debiased squared weighted phase lag index; 
%       requires the complex cross spectrum with dimensions
%       (nchan x nchan x nepochs)
%
%   - SMILe_PPC: returns the unweighted and weighted forms of the pairwise 
%       phase consistency; requires the complex cross spectrum with 
%       dimensions (nchan x nchan x nepochs)
%
%--------------------------------------------------------------------------
% Thresholding and - SMILe_CIcalc
%
% Transformation - SMILe_Fisher_Z
%--------------------------------------------------------------------------
% Each of the above functions (except for SMILe_instantPLV) returns a
% single symmetrical, real-valued, "weighted" adjacency matrix that 
% represents all possible pairwise connections between channels. 
%
% SMILe_CIcalc can be used to calculate a 95% confidence interval, which is
% returned along with standard deviation values (mean +/- sd*1, 
% mean +/- sd*2 and mean +/- sd*3). Both of these can be used as needed,
% but will principlally be used for thresholding plots.
%
% Similarly, it is often a good idea to transform connectivity data due to
% non-normal distributions and a bias in statistical calculations.
% SMILe_Fisher_Z takes an adjacency matrix (weighted) and computes the
% Fisher Z-transform on the connectivity values to approximate a normal
% distribution and remove bias.
%
%--------------------------------------------------------------------------
% Permutation test against random connections - SMILe_Perm_Zthresh_Adj
%
% Permutation test for two groups/conditions - SMILe_Perm_Zthresh_Diff_Adj
%
%--------------------------------------------------------------------------
% Plotting Adjacency Matrices
%--------------------------------------------------------------------------
% There are a few plotting options available to display connectivity in
% various ways. Each depends on an adjacency matrix (from above) and a
% channel location variable to work.
%
%   - SMILe_Line_Plot: More or less the standard connectivity plot, where
%       connecting lines are drawn between electrodes whose connection
%       strength exceeds a threshold value (or all connections if not
%       thresholded). 
%   
%       Requires an adjacency matrix (weighted or binary) and layout 
%       (from 64layout.mat). Can use SMILe_Line_Plot([],[]) for a demo.
%
%   - SMILe_Line_Plot_3D: Similar to the SMILe_Line_Plot, but results in a
%       3-D representation of electrodes and their connectivity pattern. 
%       
%
%       Requires an adjacency matrix (weighted or binary) and elec 
%       (from 64elec.mat).
%
%   - SMILe_Head_In_Head_Plot: Basically copied from Giudo Nolte and Arne
%       Ewald. Shows connectivity across the scalp, where each electrode is
%       represented and all regional connectivities (above a threshold, or 
%       not) are shown as a colored topography. The relative location of
%       the individual topographies directly correcponds with the location
%       of each electrode as it sits on the head. 
% 
%       Requires an adjacency matrix (weighted), chanlocs_2D (from 
%       chanlocs_2D_3D.mat) and a pars structure for additional settings 
%       (pars can be made empty ([]) and passed in to set default values).
%
%   - SMILe_Circle_Plot: Similar to plots that might be seen in large-scale
%       network analyses, each electrode name is displayed around the 
%       circumference of the circle, and lines are drawn between electrodes 
%       whose connection strength exceeds a threshold value (or all 
%       connections if not thresholded). 
%
%       Requires an adjacency matrix (weighted or binary) and chanlocs 
%       (from 64chanlocs.mat). Can use SMILe_Line_Plot([],[]) for a demo.
%
%--------------------------------------------------------------------------
% Graph Analysis - still under development
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------





%--------------------------------------------------------------------------
% Example usage:
% clear cs nave cfg
%
% cfg = [];
% cfg.segleng = [];
% cfg.segshiftval = [];
% cfg.segstart = [];
% cfg.fs = [];
% cfg.wintype = [];
% cfg.zeropad = [];
% cfg.subave = [1];
% cfg.segave = [1];
% cfg.fftdetrend = [];
% cfg.keeptrials = [];
% cfg.maxfreq = [];
% cfg.band = [8 12];
% 
% [cs, fftdataout] = SMILe_CS_Fourier(data, cfg);