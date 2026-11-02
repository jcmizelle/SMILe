function FAA = SMILe_FAA(data, cfg)
%%
%--------------------------------------------------------------------------
% Example user input:
%--------------------------------------------------------------------------
%cfg = [];
%cfg.Lchans = [3 4 5];
%cfg.Lchans = {'C1'; 'C3'};
%cfg.Rchans = [10 11 12];
%cfg.Rchans = {'C2'; 'C4'};
%cfg.frex = [8 12];
%cfg.srate = 1000;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Get channel location
%--------------------------------------------------------------------------
load('path/to/templates/64layout.mat')
load('path/to/templates/64chanlocs.mat')

%--------------------------------------------------------------------------
% Get label information
%--------------------------------------------------------------------------
labels = {chanlocs.labels}.';

%--------------------------------------------------------------------------
% Preallocate variables:
%--------------------------------------------------------------------------
Lchans = [];
Rchans = [];
Ldata = [];
Rdata = [];

%--------------------------------------------------------------------------
% Pull channel information (whether already indices or names - handle it
% here:
%--------------------------------------------------------------------------
if ~iscell(cfg.Lchans)
    Ldata = data(cfg.Lchans,:,:);
    Rdata = data(cfg.Rchans,:,:);
   
elseif iscell(cfg.Lchans)
    for nchans = 1:numel(cfg.Lchans)
        Lchans(nchans) = find(strcmp(labels, cfg.Lchans(nchans)));
        Rchans(nchans) = find(strcmp(labels, cfg.Rchans(nchans)));
        
    end

    Ldata = data(Lchans,:,:);
    Rdata = data(Rchans,:,:);
    
end

%--------------------------------------------------------------------------
% Handle data:
%--------------------------------------------------------------------------
dims = size(Ldata);
if numel(dims) > 2
    Ldata = mean(Ldata, 3);
    Rdata = mean(Rdata, 3);
    
end

if dims(1) > 1
    Ldata = mean(Ldata, 1);
    Rdata = mean(Rdata, 1);
    
end

%--------------------------------------------------------------------------
% Identify frequency range:
%--------------------------------------------------------------------------
FREQ_1 = cfg.frex(1);
FREQ_2 = cfg.frex(2);

%--------------------------------------------------------------------------
% Compute power spectrum for left channel(s):
%--------------------------------------------------------------------------
WIND = hamming(floor(length(Ldata))/1.5);   % Get 1.5 sec time windows
OVER = floor((length(Ldata))/1.5/2);        % 50% overlap
SIGN = Ldata';                              % Get signal
[~, freqs, ~, power] = spectrogram(SIGN, WIND, OVER, [], cfg.srate);
indFreqs = find(freqs>FREQ_1 & freqs<FREQ_2);
POW_L = power(indFreqs);

%--------------------------------------------------------------------------
% Compute power spectrum for right channel(s):
%--------------------------------------------------------------------------
WIND = hamming(floor(length(Rdata))/1.5);   % Get 0.5 sec time windows
OVER = floor((length(Rdata))/1.5/2);        % 50% overlap
SIGN = Rdata';                              % Get signal
[~, freqs, ~, power] = spectrogram(SIGN, WIND, OVER, [], cfg.srate);
indFreqs = find(freqs>FREQ_1 & freqs<FREQ_2);
POW_R = power(indFreqs);

%--------------------------------------------------------------------------
% Compute whole FAA:
%--------------------------------------------------------------------------
FAA = mean(abs(log(POW_R)-log(POW_L)));

%--------------------------------------------------------------------------
% End for function:
%--------------------------------------------------------------------------
end