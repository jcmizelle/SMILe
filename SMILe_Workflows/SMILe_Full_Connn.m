%% SMILe_Full_Conn
% 
% This script brings in preprocessed data (SMILe_Proprocessing) which
% should be clean data in the format of nchan x ntimes x ntrials
%
% This has been built to use CSD data (Laplacian transformed) and will
% check to see if this is present in the EEG structure. If not, this will
% be computed and then moved forward through the rest of the script.
%
% This function will locally execute the Laplacian (if needed), and then
% compute functional connectivity, graph theory and directed phase transfer
% entropy.
%
% Input:
% Preprocessed data from SMILe_Proprocessing in the format of 
% nchan x ntimes x ntrials
%
% Output:
% Thresholded functional connectivity
% All of the graph measures found in SMILe_GraphMeasures
% Thresholded directed phase transfer entropy
%
%%
%--------------------------------------------------------------------------
% We will eventually need electrode templates - load them here:
%--------------------------------------------------------------------------
load('path/to/templates/64elec.mat')
load('path/to/templates/64layout.mat')

%% 
%--------------------------------------------------------------------------
% Begin to identify data and set initial parameters
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Identify the data folder:
%--------------------------------------------------------------------------
dataFolder = 'path/to/data';

%--------------------------------------------------------------------------
% Identify a specific file, or compile a list depending on the search 
% parameters:
%--------------------------------------------------------------------------
files = [];
files = dir(fullfile(dataFolder,'Y_S002-EEG_2*.set'));

%--------------------------------------------------------------------------
% Sampling frequency
%--------------------------------------------------------------------------
fs = 1000;

%--------------------------------------------------------------------------
% Specific bandwidth of interest
%--------------------------------------------------------------------------
band = [8:12];
cs = [];
tmpcs = [];

%--------------------------------------------------------------------------
% Identify a time range for the analysis
%--------------------------------------------------------------------------
twindow = [1650:2250];

%% 
%--------------------------------------------------------------------------
% Bring in selected data, compute CSD (Laplacian) if not pre-existing
%--------------------------------------------------------------------------
fpath = files(nsub).folder;
fname = files(nsub).name(1:end-4);
fext = files(nsub).name(end-3:end);
EEG = pop_loadset('filename', strcat(fname, fext), 'filepath', fpath);

%--------------------------------------------------------------------------
% Check for CSD data, run CSD (Laplacian) if not present
%--------------------------------------------------------------------------
try 
    data = EEG.CSD;

catch
%--------------------------------------------------------------------------
% Grab the search path and find instances of the Laplacian code folder:
%--------------------------------------------------------------------------
    pathCell = regexp(path, pathsep, 'split');
    for p = 1:size(pathCell, 2)
      onPath(p) = any(strfind(pathCell{p}, 'SMILe_Laplace'));

    end

%--------------------------------------------------------------------------
% Make sure that the Laplacian code is fully installed:
%--------------------------------------------------------------------------
    if sum(onPath) > 1
        [CSD, ~] = CSD_Analysis(EEG);
        data = CSD;
%--------------------------------------------------------------------------
% If not, install it to the search path:
%--------------------------------------------------------------------------
    else 
        addpath(genpath(cfg.LapPath));
        [CSD, ~] = CSD_Analysis(EEG);
        data = CSD;
        
    end

end

%%
%--------------------------------------------------------------------------
% Perform a bandpass filter for the frequency range specified above.
%
% The SMILe_Wavelet_BandFilter function will export several things, but
% here we are spceifically pulling the filteres time series as well as the
% crossspectrum
% 
% This also collapses the crossspectrum within a time window, as set above,
% but leaves the filtered time series intact.
%--------------------------------------------------------------------------
[fdata, ~, ~, cs] = SMILe_Wavelet_BandFilter(data,fs,band,1);

tmpcs(:,:,:) = squeeze(mean(cs(:,:,twindow,:),3));

%% 
%--------------------------------------------------------------------------
% Compute connectivity - here with dswpli, plenty of other options
%--------------------------------------------------------------------------
adj = SMILe_wPLI(tmpcs);
adj = triu(adj);
adj = adj + adj';
%adj(adj~=0) = 1;
adj(logical(eye(size(adj)))) = 0;
adj = abs(adj);
%adj(adj == 1) = .999;
%adj = atanh(adj);

%--------------------------------------------------------------------------
% Threshold for adjacency matrices
%--------------------------------------------------------------------------
cfg = [];
cfg.nperms = 1000;
cfg.pval = 0.05;
cfg.doplot = 1;
cfg.layout = layout;

[thr_adj] = SMILe_Perm_Zthresh_Adj(adj, cfg);

%% 
%--------------------------------------------------------------------------
% Calculate MST and Graph Measures
%
% ECOfilter is used here. Lots of other options.
%--------------------------------------------------------------------------
mst_adj = zeros(size(adj));
tmp_adj = SMILe_ECOfilter(adj,0);
nz = find(tmp_adj ~= 0);
mst_adj(nz) = adj(nz);
mst_adj = tmp_adj;
mst_adj = triu(mst_adj);
mst_adj = mst_adj + mst_adj';
%mst_adj(adj~=0) = 1;
mst_adj(logical(eye(size(mst_adj)))) = 0;

%% 
%--------------------------------------------------------------------------
% Construct the graph for visualization in a circuit diagram
%--------------------------------------------------------------------------
[G] = SMILe_adj2graph_binary(adj,elec);
SMILe_Circuit_Plot(G, 1);
SMILe_Circuit_Plot_MST(mst_adj, elec);

%% 
%--------------------------------------------------------------------------
% Extract the graph measures
%--------------------------------------------------------------------------
GrMeas = SMILe_GraphMeasures(mst_adj);
GrMeas.Ecc(isnan(GrMeas.Ecc)) = 0;

pars = [];
%pars.scale = [];
pars.cmap = 'jet';
pars.resolution = 300;
pars.cbar = 1;
SMILe_TopoPlot(GrMeas.Ecc, layout, pars);
SMILe_TopoPlot_Hubs(layout, pars, GrMeas.hubs)

%% 
%--------------------------------------------------------------------------
% Compute the directed phase transfer entropy:
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Set another time window if desired (should leave the same as above,
% though)
%--------------------------------------------------------------------------
%twindow = [1700:2000];

%--------------------------------------------------------------------------
% Do the calculations:
%--------------------------------------------------------------------------
[dPTE, ~] = SMILe_PhaseTE(fdata(:,twindow,:), [], 'scott');

%--------------------------------------------------------------------------
% Stastically threshold based on surrogate data. Here we are using the AAFT 
% process. Lots of other options.
%
% Basically repeat the whole dPTE process after computing the surrogates.
%
% Loop to compute dPTE on each surrogate repetition for thresholding.
%--------------------------------------------------------------------------
%[surr] = SMILe_Surrogates_RandPhase(data(:,twindow,:), 100);
[surr] = SMILe_Surrogates_AAFT(data(:,twindow,:), 1000);
%[surr] = SMILe_Surrogates_iAAFT(data(:,twindow,:), 100, 50);

[sfdata, ~, ~, ~] = SMILe_Wavelet_BandFilter(surr,fs,band,0);

rpt = 'Surrogate dPTE';
for ns = 1:size(sfdata,3)
    SMILe_TextProgressBar(ns, size(sfdata,3), rpt);
    [outdPTE, ~] = SMILe_PhaseTE(sfdata(:,:,ns), [], 'scott');
    surrogdPTE(:,:,ns) = outdPTE;
    clear outdPTE

end

%--------------------------------------------------------------------------
% Actually apply the threshold here:
%--------------------------------------------------------------------------
[thrdPTE] = SMILe_Surrogates_Threshold(dPTE, surrogdPTE, 1, 0.1);

%--------------------------------------------------------------------------
% Plot if desired
%--------------------------------------------------------------------------
SMILe_Line_Arrow_Plot(thrdPTE, layout, 0, 1);

%--------------------------------------------------------------------------
% End of file
%--------------------------------------------------------------------------
