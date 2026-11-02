%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ SMILe Preprocessing Routine - Configuration Details
%$$$$$$ JCM, current as of 05/11/2021
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$ Required Settings
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
cfg = [];

%--------------------------------------------------------------------------
% Do you want to re-process data? reRun = 0 (No) or reRun = 1 (Yes). If
% reRun is set to 0, already processed files will be skipped.
%--------------------------------------------------------------------------
cfg.reRun = 0;

%--------------------------------------------------------------------------
% Data folder - the actual location of the data.
%--------------------------------------------------------------------------
cfg.dataFolder = 'path/to/data';

%--------------------------------------------------------------------------
% Do you want to process all files in a folder, or only select a subset of 
% the files? Set searchString = []; to process all files. Set 
% searchString = 'S001'; to process only the siles associated with this
% participant/condition/etc.
%--------------------------------------------------------------------------
% cfg.searchString = [];
cfg.searchString = 'O*'; 

%--------------------------------------------------------------------------
% What type of files are you working with?
%--------------------------------------------------------------------------
cfg.fileType = '.cnt';

%--------------------------------------------------------------------------
% Where should output go when complete?
%--------------------------------------------------------------------------
cfg.outFolder = 'path/to/data';

%--------------------------------------------------------------------------
% The Laplacian transform is always calculated. It is not a plugin, so the
% location of the base code should be provided.
%--------------------------------------------------------------------------
cfg.LapPath = 'path/to/SMILe_Laplace';

%--------------------------------------------------------------------------
% Are there multiple conditions to consider in the file?
%--------------------------------------------------------------------------
cfg.multconds = 1;

%--------------------------------------------------------------------------
% If so, what are the event codes embedded in the file?
%--------------------------------------------------------------------------
cfg.evmarkers = [2 4 8 16 32 64];

%--------------------------------------------------------------------------
% Would you like for the unique conditions to have specific output names, 
% instead of using the event codes? Note: the order of the optional names
% must match the order of the event codes!
%--------------------------------------------------------------------------
cfg.conds = {'Cor', 'Plaus', 'Cont', 'Manip', 'Hand', 'Tool'};

%--------------------------------------------------------------------------
% Are there events in the data set that you know you want to remove (if you
% do not specify event codes above)?
%--------------------------------------------------------------------------
cfg.killevmarkers = 255;

%--------------------------------------------------------------------------
% Is this a continuous data set with no events?
%--------------------------------------------------------------------------
cfg.continuous = 0;

%--------------------------------------------------------------------------
% Optionally add pesudo events into continuous data.
% addEvents = 0 (No) or addEvents = 1 (Yes). Default is no.
%--------------------------------------------------------------------------
cfg.addEvents = 0;

%--------------------------------------------------------------------------
% If the above is true, events will automatically be inserted in the file. 
% The default epoch length is 2, but can be changed here:
%--------------------------------------------------------------------------
cfg.contEpLeng = 2;

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$ Optional Settings
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

%--------------------------------------------------------------------------
% Do you want to change the sampling rate of the data? newSampleRate = [];
% will not alter the sampling rate, but newSampleRate = 128; will resample
% to 128Hz, for example.
%--------------------------------------------------------------------------
cfg.newSampleRate = [];

%--------------------------------------------------------------------------
% Do you want to remove channels from the montage? Identify their indices
% here. killchanids = []; for no channel removal.
%--------------------------------------------------------------------------
cfg.killchanids = [67 68 69];

%--------------------------------------------------------------------------
% Do you want to remove channels from the montage but not permanently lose 
% them? Identify their indices names here. movechanids = []; for no channel 
% removal.
%--------------------------------------------------------------------------
cfg.movechanids = {'M1', 'M2', 'VEOG', 'HEOG'};

%--------------------------------------------------------------------------
% High pass filter cutoff.
%--------------------------------------------------------------------------
cfg.hpFilterCut = 1;

%--------------------------------------------------------------------------
% Maximum frequency of interest (for removal of line noise if necessary).
%--------------------------------------------------------------------------
cfg.maxfreq = 40;

%--------------------------------------------------------------------------
% Low pass filter cutoff.
%--------------------------------------------------------------------------
cfg.lpFilterCut = 35;

%--------------------------------------------------------------------------
% Do you want to clean your data with Artifact Subspace Reconstruction? 
% doASR = 0 (No) or doASR = 1 (Yes). Default is yes.
%--------------------------------------------------------------------------
cfg.doASR = 1;

%--------------------------------------------------------------------------
% Do you want re-reference using an infinity reference? 
% doREST = 0 (No) or doREST = 1 (Yes). Default is no - the average 
% reference is computed by default.
%--------------------------------------------------------------------------
cfg.doREST = 0;

%--------------------------------------------------------------------------
% REST is not currently a plugin - need to identify the location for the 
% base code.
%--------------------------------------------------------------------------
cfg.RESTpath = 'path/to/SMILe_Infinity_Reference';

%--------------------------------------------------------------------------
% Optionally return previously moved channels to the montage.
% moveBack = 0 (No) or moveBack = 1 (Yes). Default is no.
%--------------------------------------------------------------------------
cfg.moveBack = 0;

%--------------------------------------------------------------------------
% For event-related data. Start of an epoch (in seconds). Negative time
% reflects pre-stimulus time.
%--------------------------------------------------------------------------
cfg.epstart = -3;

%--------------------------------------------------------------------------
% For event-related data. End of an epoch (in seconds). Positive time
% reflects post-stimulus time.
%--------------------------------------------------------------------------
cfg.epend = 3;

%--------------------------------------------------------------------------
% Perform Independent Components Analysis?
% doICA = 0 (No) or doICA = 1 (Yes). Default is no.
%--------------------------------------------------------------------------
cfg.doICA = 1;

%--------------------------------------------------------------------------
% If ICA is to be calculated, what method?
% AMICA: ICAtype = 1; 
% Infomax: ICAtype = 2; 
% SOBI: ICAtype = 3; 
% JADE: ICAtype = 4; 
% PICARD: ICAtype = 5; (Default)
%--------------------------------------------------------------------------
cfg.ICAtype = 5;

%--------------------------------------------------------------------------
% Remove components identified as artiface?
% removeComponents = 0 (No) or removeComponents = 1 (Yes). Default is no.
%--------------------------------------------------------------------------
cfg.removeComponents = 0;

%--------------------------------------------------------------------------
% Identify components whether or not you remove them?
% classifyComponents = 0 (No) or removeComponents = 1 (Yes). Default is 
% yes.
%--------------------------------------------------------------------------
cfg.classifyComponents = 1;

%--------------------------------------------------------------------------
% If component removal is selected, what classifier? 
% ICLabel: CLASStype = 1;
% SASICA: CLASStype = 2;
% MARA: CLASStype = 3;
% ADJUST: CLASStype = 4;
%--------------------------------------------------------------------------
cfg.CLASStype = 1;

%--------------------------------------------------------------------------
% For event-related data. Potentially new start of an epoch (in seconds). 
% Negative time reflects pre-stimulus time.
%--------------------------------------------------------------------------
cfg.newepstart = -1.5;

%--------------------------------------------------------------------------
% For event-related data. Potentially new end of an epoch (in seconds). 
% Positive time reflects post-stimulus time.
%--------------------------------------------------------------------------
cfg.newepend = 2;

%--------------------------------------------------------------------------
% For event-related data. Baseline removal for epochs (in milliseconds). 
% Time starts at the beginning of an epoch, which can be negative time, but 
% this value needs to be positive..
%--------------------------------------------------------------------------
cfg.baselength = 500;

%--------------------------------------------------------------------------
% Optionally compute dipoles if ICA was calculated.
% fitDipoles = 0 (No) or fitDipoles = 1 (Yes). Default is no.
%--------------------------------------------------------------------------
cfg.fitDipoles = 0;

%--------------------------------------------------------------------------
% If dipoles are computed, which method?
% estimationType = 1; for BESA estimationType = 2; for BEM.
%--------------------------------------------------------------------------
cfg.estimationType = 2;

%--------------------------------------------------------------------------
% Optionally compute time-frequency properties.
% calcTF = 0 (No) or calcTF = 1 (Yes). Default is no.
%--------------------------------------------------------------------------
cfg.calcTF = 0;

%--------------------------------------------------------------------------
% Optionally compute spectral properties.
% doSpec = 0 (No) or doSpec = 1 (Yes). Default is no.
%--------------------------------------------------------------------------
cfg.doSpec = 0;

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$ Execute
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
SMILe_Preprocessing(cfg)