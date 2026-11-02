function SMILe_Video_Maker(datain, cfg)
%--------------------------------------------------------------------------
%
% This function produces a topographical plot of activation (any type)
% across the electrode montage.
%
% This is an improved version of the EEGLAB topoplot, and mimics the
% FieldTrip version without all of the headache.
%
% Input:
%   datain      - A nchan x 1 vector of data to visualize
%
%   cfg         - A configuration variable
%   cfg.datatimes:      Data points of interest (in samples)
%   cfg.timeStamp:      Add a time stamp. 0 or 1. Default is 1.
%   cfg.realtimes:      Real times (in ms) - for including a time stamp
%   cfg.timezero:       Time of stimulus presentation (in samples)
%   cfg.timestep:       Time resolution
%   cfg.workingDir:     Location for images to be saved
%   cfg.outdir:         Location for the video to be saved
%   cfg.videoName:      Name for the video to be saved
%   cfg.frameRate:      Frame rate for the video
%   cfg.quality:        Quality of the saved video
%
%   cfg.pars     - A configuration variable, where: 
%   pars.scale:         Sets the scale of the color map. Either a 1x2 
%                       vector corresponding to minimum and a maximum or 
%                       just a number (say x) then the scale is from 
%                       [-x x]. The default is scale = [min(z),max(z)];
%   pars.resolution:    Sets the resolution - default is 300
%   pars.cbar:          Add a colorbar. 1 or 0; Default it 0
%   pars.cbarsize:      Scales the size of the colorbar by this factor.
%                       Default is .9;
%   pars.addContour:    Add contour lines to the plot. 1 or 0. Default is
%                       0.
%   pars.addElec:       Add electrodes to the plot; 0 or 1. Default is 1.
%
%--------------------------------------------------------------------------
%%
%--------------------------------------------------------------------------
% Pull configurations from the cfg:
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Times of interest
%--------------------------------------------------------------------------
datatimes = cfg.datatimes;

%--------------------------------------------------------------------------
% Times of interest
%--------------------------------------------------------------------------
realtimes = cfg.realtimes;

%--------------------------------------------------------------------------
% Time resolution
%--------------------------------------------------------------------------
timestep = cfg.timestep;

%--------------------------------------------------------------------------
% Time of stimulus
%--------------------------------------------------------------------------
timezero = cfg.timezero;

%--------------------------------------------------------------------------
% Output for the images
%--------------------------------------------------------------------------
workingDir = cfg.workingDir;

%--------------------------------------------------------------------------
% Directory for the video to be saved
%--------------------------------------------------------------------------
outdir = cfg.outdir;

%--------------------------------------------------------------------------
% Name for the video to be saved under
%--------------------------------------------------------------------------
videoName = cfg.videoName;

%--------------------------------------------------------------------------
% Frame rate of the video
%--------------------------------------------------------------------------
try 
    frameRate = cfg.frameRate;
    
catch
    frameRate = 15;
    
end

%--------------------------------------------------------------------------
% Video resolution
%--------------------------------------------------------------------------
try
    quality = cfg.quality;
    
catch
    quality = 100;
    
end

%--------------------------------------------------------------------------
% Parameters for the topo plot
%--------------------------------------------------------------------------
try
    pars = cfg.pars;
    
catch
    pars = [];
    
end

%--------------------------------------------------------------------------
% Layout for the topo plot:
%--------------------------------------------------------------------------
try
    layout = cfg.layout;
    
catch
    load('path/to/templates/64layout')
    
end

%%
%--------------------------------------------------------------------------
% Check the existence of the output folder, and create it if needed
% Change current directory to the output folder
%--------------------------------------------------------------------------
cwd = pwd;

if ~exist('workingDir','dir') == 1
    mkdir(char(workingDir))
   
end

cd(workingDir)

%--------------------------------------------------------------------------
% Build the timeline for "temporal resolution" - one image will be created
% for each time step
%--------------------------------------------------------------------------
tline = [datatimes(1):timestep:datatimes(2)];

%--------------------------------------------------------------------------
% Build the timeline for "temporal resolution" - one image will be created
% for each time step
%--------------------------------------------------------------------------
realtime = [realtimes(1):timestep:realtimes(2)];

%--------------------------------------------------------------------------
% Prepare the data (average with time steps) and plot
% Also save the images to the desired outpot folder
%--------------------------------------------------------------------------
pars.scale = [min(min(datain))*.60 max(max(datain))*.60];
%pars.addContour = 1;
for i = 2:size(tline,2)
    pdata = mean(datain(:,[tline(i-1):tline(i)]), 2);
    SMILe_TopoPlot(pdata, layout, pars)
    if cfg.timeStamp == 1
        SMILe_figtitle([num2str(realtime(i)) 'ms'],'fontweight','bold','fontsize',30)
    
    end
    
    saveas(gcf,sprintf('FIG-SVM-%03d.jpg',i))
    close
    pause(.2)

end

%%
%--------------------------------------------------------------------------
% Begin the video making process
% Also save the images to the desired outpot folder
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Collect names of all of the images
%--------------------------------------------------------------------------
imageNames = dir(fullfile(workingDir,'FIG-SVM*.jpg'));
imageNames = {imageNames.name}';

%--------------------------------------------------------------------------
% Specify video frame rate and quality
% Initialize the video maker
%--------------------------------------------------------------------------
outdir = fullfile(workingDir, 'Video');

if ~exist(outdir, 'dir')
    mkdir(char(outdir))
    
end

outputVideo = VideoWriter(fullfile(workingDir, 'Video', videoName),'MPEG-4');
outputVideo.FrameRate = frameRate;
optputVideo.Quality = quality;
open(outputVideo)

%--------------------------------------------------------------------------
% Load images, concatenating into a video
%--------------------------------------------------------------------------
for ii = 1:length(imageNames)
   img = imread(fullfile(workingDir, imageNames{ii}));
   writeVideo(outputVideo, img)
   
end

%--------------------------------------------------------------------------
% Close video maker
%--------------------------------------------------------------------------
close(outputVideo)

%--------------------------------------------------------------------------
% Move pack to previous directory
%--------------------------------------------------------------------------
cd(cwd)

%--------------------------------------------------------------------------
% End of function
%--------------------------------------------------------------------------
end