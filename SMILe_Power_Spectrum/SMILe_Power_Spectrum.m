function powspect = SMILe_Power_Spectrum(EEG)
dataFolder = 'path/to/data';

files = dir(fullfile(dataFolder, '*.set'));

%outFolder = fullfile(fpath, 'CS');
outFolder = 'path/to/data';

if ~exist('outFolder','dir') == 1
    mkdir(char(outFolder))
   
end

for i = 1%1:size(files,1)
fpath = files(i).folder;

fname = files(i).name(1:end-4);

fext = files(i).name(end-3:end);

EEG = pop_loadset('filename',strcat(fname, fext),'filepath',fpath);


end

%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Assign locations for all electrodes:
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
tmpdipfitfolder = dir([eeglabdir filesep 'plugins' filesep 'dipfit*']);
dipfitfolder = fullfile(tmpdipfitfolder.folder, tmpdipfitfolder.name);
clear tmpdipfitfolder

%--------------------------------------------------------------------------
% Try for a Windows machine:
%--------------------------------------------------------------------------
try EEG = pop_chanedit(EEG, 'lookup',...
    [dipfitfolder '\standard_BESA\standard-10-5-cap385.elp']);

%--------------------------------------------------------------------------
% Catch for Linux and Unux machines:
%--------------------------------------------------------------------------
catch
    EEG = pop_chanedit(EEG, 'lookup',...
    [dipfitfolder '/standard_BESA/standard-10-5-cap385.elp']);

end

%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Run a highpass filter to detrend the data, removing very slow
%$$$$$$ potentials and amplifier drift.
%$$$$$$
%$$$$$$ Syntax: pop_eegfiltnew(EEG, locutoff, hicutoff, filtorder,revfilt, usefft, plotfreqz, minphase)
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
if isfield(cfg, 'hpFilterCut') && ~isempty(cfg.hpFilterCut)
    hpFilterCut = cfg.hpFilterCut;

else
    hpFilterCut = 1;
    
end

%--------------------------------------------------------------------------
% Do the filtering:
%--------------------------------------------------------------------------
EEG = pop_eegfiltnew(EEG, [], hpFilterCut, [], true, [], 0);

%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Run a lowpass filter to eliminate (high) frequeicies of no interest
%$$$$$$
%$$$$$$ Syntax: pop_eegfiltnew(EEG, locutoff, hicutoff, filtorder,revfilt, usefft, plotfreqz, minphase)
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
if isfield(cfg, 'lpFilterCut') && ~isempty(cfg.lpFilterCut)
    lpFilterCut = cfg.lpFilterCut;

else
    lpFilterCut = 35;
    
end

%--------------------------------------------------------------------------
% Do the filtering:
%--------------------------------------------------------------------------
EEG = pop_eegfiltnew(EEG, [], lpFilterCut, [], false, [], 0);

originalEEG = EEG;

%%
EEG = clean_artifacts(EEG,...
                            'ChannelCriterion', 0.85,...
                            'LineNoiseCriterion', 4,...
                            'BurstCriterion', 5,...
                            'WindowCriterion', 0.275,...
                            'Highpass', [],...
                            'FlatlineCriterion', 5,...
                            'BurstRejection', 'off',...
                            'Distance', 'euclidian',...
                            'MaxMem', 8192);
                        

%%
if EEG.nbchan ~= originalEEG.nbchan
    EEG = pop_interp(EEG, originalEEG.chanlocs, 'spherical');

end

%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Add path(s) for Laplacian transform.
%$$$$$$
%$$$$$$ I always include this as a final strp and include the resulting
%$$$$$$ computations in the EEG structure.
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%addpath(genpath([eeglabdir filesep 'Laplace' filesep]));
%rehash
if isfield(cfg, 'LapPath') && ~isempty(cfg.LapPath)
    LapPath = cfg.LapPath;
    
else LapPath = uigetdir;
    
end

addpath(genpath(LapPath));
rehash
%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Compute the Laplacian transformation of the data and add to the EEG
%$$$$$$ structure, but not replacing the EEG.data structure.
%$$$$$$
%$$$$$$ Both measures are based on the Perrin method and yield similar
%$$$$$$ results. The CSD function also returns the spherical spline
%$$$$$$ surface potential (SP) interpolation.
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

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
   [CSDout, SLAP] = laplacian_CSD(EEG);
   EEG.CSD = CSDout;
   EEG.SLAP = SLAP;

%--------------------------------------------------------------------------
% If not, install it to the search path:
%--------------------------------------------------------------------------
else 
   addpath(genpath(cfg.LapPath));
   [CSDout, SLAP] = laplacian_CSD(EEG);
   EEG.CSD = CSDout;
   EEG.SLAP = SLAP;
   
end

%% Get the spectral data

%--------------------------------------------------------------------------
% Run spectopo for every electrode:
%--------------------------------------------------------------------------
for channelIdx = 1:EEG.nbchan
    [psdOutDb(channelIdx,:), freq] = spectopo(EEG.CSD(channelIdx, 1:500), 0, EEG.srate, 'plot', 'off');

%--------------------------------------------------------------------------
% Get the indices of the frequency bands of interest:
%--------------------------------------------------------------------------
deltaIdx        = find(freqs>1 & freqs<4);
thetaIdx        = find(freqs>4 & freqs<8);
alphaIdx        = find(freqs>8 & freqs<12);
muIdx           = find(freqs>10 & freqs<12);
sigmaIdx        = find(freqs>12 & freqs<15);
fullbetaIdx     = find(freqs>13 & freqs<30);
lowbetaIdx      = find(freqs>13 & freqs<18);
highbetaIdx     = find(freqs>18 & freqs<30);
smbetaIdx       = find(freqs>18 & freqs<22);

%--------------------------------------------------------------------------
% Compute absolute power:
%--------------------------------------------------------------------------
deltaPower(channelIdx)      = mean(10.^(psdOutDb(deltaIdx)/10));
thetaPower(channelIdx)      = mean(10.^(psdOutDb(thetaIdx)/10));
alphaPower(channelIdx)      = mean(10.^(psdOutDb(alphaIdx)/10));
muPower(channelIdx)         = mean(10.^(psdOutDb(muIdx)/10));
sigmaPower(channelIdx)      = mean(10.^(psdOutDb(sigmaIdx)/10));
fullbetaPower(channelIdx)   = mean(10.^(psdOutDb(fullbetaIdx)/10));
lowbetaPower(channelIdx)    = mean(10.^(psdOutDb(lowbetaIdx)/10));
highbetaPower(channelIdx)   = mean(10.^(psdOutDb(highbetaIdx)/10));
smbetaPower(channelIdx)     = mean(10.^(psdOutDb(smbetaIdx)/10));

%--------------------------------------------------------------------------
% Convert to dB:
%--------------------------------------------------------------------------
deltadB(channelIdx)      = 10*log10(deltaPower(channelIdx));
thetadB(channelIdx)      = 10*log10(thetaPower(channelIdx));
alphadB(channelIdx)      = 10*log10(alphaPower(channelIdx));
mudB(channelIdx)         = 10*log10(muPower(channelIdx));
sigmadB(channelIdx)      = 10*log10(sigmaPower(channelIdx));
fullbetadB(channelIdx)   = 10*log10(fullbetaPower(channelIdx));
lowbetadB(channelIdx)    = 10*log10(lowbetaPower(channelIdx));
highbetadB(channelIdx)   = 10*log10(highbetaPower(channelIdx));
smbetadB(channelIdx)     = 10*log10(smbetaPower(channelIdx));

clear psdOutDb
%---------------------------------------------------------------------
% End loop for electrodes:
%--------------------------------------------------------------------------
end

%% EGGHead Plot (stolen from ERPLAB)
addpath('path/to/SMILe_EggHeadPlot');
addpath('path/to/SMILe_Plotting');

channelvector = {layout.label};

topodata = fullbetadB;

SMILe_eggheadplot(...
    'Channels', channelvector{:}, ...
    'Amplitude', topodata, ...
    'Knockout', [], ...
    'Method', 'v4', ...
    'Scale', [min(topodata) max(topodata)], ...
    'Contours', 0, ...
    'MapStyle', 'jet', ...
    'Steps', [512], ...
    'Smooth', [0], ...    
    'Style', 'Full', ...
    'ShowElectrodes', 'Yes', ...
    'ElectrodeSize', 20, ...
    'FontSize', 10, ...
    'Label', [], ...
    'ShowBrain', 'Yes', ...
    'BrainOpacity', 0.5 ...
    );


%% My topoplot
%--------------------------------------------------------------------------
% Load a channel configuration files:
%--------------------------------------------------------------------------
load('path/to/templates/32elec.mat');
load('path/to/templates/32layout.mat');

%--------------------------------------------------------------------------
% Define topo data:
%--------------------------------------------------------------------------
topodata = thetadB;
topodata = alphadB;
topodata = fullbetadB;
topodata = mudB;

%--------------------------------------------------------------------------
% Create the plot:
%--------------------------------------------------------------------------
pars = [];
pars.scale = [min(topodata) max(topodata)];
pars.resolution = 300;
pars.cbar = [];
SMILe_TopoPlot(topodata, layout, pars)

%%
end