%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ SMILe Preprocessing Routine - Configuration Parameters
%$$$$$$ JCM, current as of 03/07/2022
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$ Required Settings
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
cfg = [];
cfg.reRun = 0;
cfg.dataFolder = 'path/to/data';
cfg.searchString = ['EXER_15']; %cfg.searchString = 'S036';
cfg.fileType = '.cdt';
cfg.doRecursive = 0;
cfg.outFolder = 'path/to/data';
cfg.doLap = 1;
cfg.LapPath = 'path/to/SMILe_Laplace';
cfg.multconds = 1;
cfg.repconds = 1;
cfg.evmarkers = [33 34];
cfg.conds = {'EXE', 'NON'};
cfg.killevmarkers = 255;
cfg.continuous = 0;
cfg.addEvents = 0;
cfg.contEpLeng = 2;

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$ Optional Settings
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
cfg.newSampleRate = [];
cfg.killchanids = [67 68 69];
cfg.movechanids = {'M1', 'M2', 'VEOG', 'HEOG'};
cfg.hpFilterCut = .1;
cfg.maxfreq = 45;
cfg.lpFilterCut = 45;
cfg.doASR = 1;
cfg.doREST = 0;
cfg.RESTpath = 'path/to/SMILe_Infinity_Reference';
cfg.moveBack = 0;
cfg.epstart = -3;
cfg.epend = 3;
cfg.doICA = 0;
cfg.ICAtype = 5;
cfg.removeComponents = 0;
cfg.classifyComponents = 0;
cfg.CLASStype = 1;
cfg.newepstart = -1;
cfg.newepend = 2;
cfg.baselength = 500;
cfg.fitDipoles = 0;
cfg.estimationType = 2;
cfg.calcTF = 0;
cfg.doSpec = 0;

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$$$$$ Execute
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
SMILe_Preprocessing(cfg)