%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Routine for computing inverse solutions from pre-processed data, as
%$$$$$$ run through the SMILe pipeline.
%$$$$$$
%$$$$$$ Here, we are using generic anatomy, since we do not have
%$$$$$$ individual MRI scans
%$$$$$$ 
%$$$$$$ This script is dependent on the EEGLAB plugin "head model", and any
%$$$$$$ and any of its dependencies, which include openmeeg
%$$$$$$
%$$$$$$ JCM, current as of 12/08/2020
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Start file
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Set path to EEGLAB parent directory
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
str = which('eeglab');
[fpath fname ext] = fileparts(str);
eeglabdir = fpath; 

%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Choose a template model for inverse solution:
%$$$$$$
%$$$$$$ Atlas_Desikan-Killiany [488]_Colin27_2003.mat                  
%$$$$$$ Atlas_Desikan-Killiany [488]_Colin27_5003.mat                  
%$$$$$$ head_modelColin27_10003_Standard-10-5-Cap339-Destrieux148.mat  
%$$$$$$ head_modelColin27_2003_Biosemi_Standard-Cap256.mat             
%$$$$$$ head_modelColin27_2003_Standard-10-5-Cap339.mat                
%$$$$$$ head_modelColin27_2003_Standard-10-5-Cap65.mat                 
%$$$$$$ head_modelColin27_5003_Standard-10-5-Cap339-Destrieux148.mat   
%$$$$$$ head_modelColin27_5003_Standard-10-5-Cap339-ROI488.mat         
%$$$$$$ head_modelColin27_5003_Standard-10-5-Cap339.mat                
%$$$$$$ head_modelColin27_8003_Standard-10-5-Cap339.mat                
%$$$$$$ head_modelColin27_8003_Standard-10-5-Cap60.mat   
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
hmodel = 'head_modelColin27_8003_Standard-10-5-Cap339.mat';
template = [eeglabdir filesep 'plugins' filesep 'dsi' filesep 'headModel' filesep 'resources' filesep hmodel];

%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ A catch to rename CB if needed
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
for loc = 1:size(EEG.chanlocs,2)
    if strcmpi(EEG.chanlocs(loc).labels, 'CB1') == 1 || strcmpi(EEG.chanlocs(loc).labels, 'I1') == 1
        EEG.chanlocs(loc).labels = 'OI1';
        dolookup = 1;

    elseif strcmpi(EEG.chanlocs(loc).labels, 'CB2') == 1  || strcmpi(EEG.chanlocs(loc).labels, 'I2') == 1
        EEG.chanlocs(loc).labels = 'OI2';
        dolookup = 1;

    elseif strcmpi(EEG.chanlocs(loc).labels, 'HEO') == 1
        EEG.chanlocs(loc).labels = 'HEOG';
        dolookup = 1;
       
    elseif strcmpi(EEG.chanlocs(loc).labels, 'VEO') == 1
        EEG.chanlocs(loc).labels = 'VEOG';
        dolookup = 1;
   
    end

end

tmpdipfitfolder = dir([eeglabdir filesep 'plugins' filesep 'dipfit*']);
dipfitfolder = fullfile(tmpdipfitfolder.folder, tmpdipfitfolder.name);
clear tmpdipfitfolder

EEG = pop_chanedit(EEG, 'lookup', [dipfitfolder filesep 'standard_BESA' filesep 'standard-10-5-cap385.elp']);


chanlocs = adjustlocs(EEG.chanlocs, 'autocenter', 'on', 'autorotate', 'on', 'autoscale', 'on');
EEG.chanlocs = chanlocs;
eeglab redraw

%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Set conductivities for the scalp, skull and brain, respectively
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
conductivity = [0.33 0.022 0.33]; 

%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Set whether to calculate the orientation free lead field (true) or 
%$$$$$$ or constrain the dipoles to be normal to the cortical 
%$$$$$$ surface (false)
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
orientation = false;              

%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Build a new head model and compute lead field matrix
%$$$$$$
%$$$$$$ This process is skipped for templates, only required for individual
%$$$$$$ anatomical scans
%$$$$$$
%$$$$$$ In either case, the resulting EEG structure will have the file 
%$$$$$$ path to the head model created in EEG.etc.src.hmfile
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
EEG = pop_forwardModel(EEG, template, conductivity, orientation, 1, 0);

%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Load and plot the computed head model, if desired:
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
hm = headModel.loadFromFile(EEG.etc.src.hmfile);
%fig = hm.plot();

%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Set parameters to calculate the inverse solution
%$$$$$$
%$$$$$$ windowSize: the number of consecutive samples used to estimate a
%$$$$$$ chunk of source time series
%$$$$$$
%$$$$$$ saveFull: whether to return the estimated time series of each
%$$$$$$ cortical dipole (true) or just the ROI reduced time series (false)
%$$$$$$
%$$$$$$ solverType: set the inverse solution method - 
%$$$$$$ 1 = loreta
%$$$$$$ 2 = Recursive Sparse Bayesian Learning
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
windowSize = 8;
saveFull = true;  
invSolver = 2;

%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Calculate the inverse solution to estimate the cortical generator
%$$$$$$ sources
%$$$$$$
%$$$$$$ The resulting EEG structure will have the sources stored in the 
%$$$$$$ EEG.etc.src field
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
switch invSolver
    case 1
        solverType = 'loreta';
        EEG = pop_inverseSolution(EEG, windowSize, saveFull, solverType);

%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Calculate the inverse solution with Recursive Sparse Bayesian 
%$$$$$$ Learning Inverse Filtering
%$$$$$$
%$$$$$$ The resulting EEG structure will have the sources stored in the 
%$$$$$$ EEG.etc.src field
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    case 2
        account4artifacts = true;
        src2roiReductionType = 'mpower';
        solverType = 'bsbl';
        updateFreq = 10;
        EEG = pop_rsbl(EEG, saveFull, account4artifacts, src2roiReductionType, solverType, updateFreq);

end
%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ End of file
%$$$$$$ 
%$$$$$$ For visualization of results, see HeadModel_LORETA_Visuals
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

