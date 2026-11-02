function [instPLV, bcinstPLV] = SMILe_instPLV(EEGdata, band, fs)
%--------------------------------------------------------------------------
%
% This function calculates time series of phase locking values, allowing 
% one to assess changes in connectivity as a function of time.
%
% Input:
% EEGdata  -    The input EEG should normally be an EEGLAB structure that
%               has been proprocessed through our lab script (epoched and
%               in the format of chan x timepoints x trials (epochs).
%
% band -        The frequency range of interest for the calculations.
%               band is numeric and must be in the range of [1 7].
%               1 - theta band (4 - 8 Hz)
%               2 - alpha band (8 - 12 Hz)
%               3 - mu band (10 - 12 Hz)
%               4 - full beta band (13 - 30 Hz)
%               5 - low beta band (13 - 18 Hz)
%               6 - high beta band (18 - 30 Hz)
%               7 - "sensorimotor" beta band (18 - 22 Hz)
%
% fs -          sampling frequency
%
% Output:
% iplv:         instaneous PLV - a time series of PLV
%
% bciplv:       Baseline Corrected instaneous PLV - a time series of
%               PLV baseline corrected per-epoch
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Declare sampling rate and band of interest:
%-------------------------------------------------------------------------- 

%--------------------------------------------------------------------------
% Redefine data. Try first to get Laplacian data, then fall back to
% non-Laplacian data. If the input is not a structure at all, just assign
% it to data.
%--------------------------------------------------------------------------
try
    if isstruct(EEGdata) && isfield(EEGdata,'CSD')
        data = EEGdata.CSD;
        disp('Using Laplacian transformed data');
    
    elseif isstruct(EEGdata)
        data = EEGdata.data;
        disp('Using non-Laplacian transformed data');
        
    end

catch
    data = EEGdata;
    disp('Data are not necessarily in EEGLAB format');
   
end

clear EEGdata

%--------------------------------------------------------------------------
% Confirm double precision:
%--------------------------------------------------------------------------
data = double(data);

%--------------------------------------------------------------------------
% Get dimensions of the data:
%--------------------------------------------------------------------------
dims = size(data);

%--------------------------------------------------------------------------
% Start progress bar:
%--------------------------------------------------------------------------
% h = waitbar(0,'Starting...');


%--------------------------------------------------------------------------
% Loop across trials:
%--------------------------------------------------------------------------
for ntrials = 1:dims(3)
   
%--------------------------------------------------------------------------
% Update progress bar:
%--------------------------------------------------------------------------
%relativeValue = nsub/size(files,1);
%waitbar(relativeValue, h, sprintf('%12.0f %%',relativeValue*100))
% waitbar(ntrials/dims(3), h, sprintf('Working on trial: %6.0f ', ntrials))
   
%--------------------------------------------------------------------------
% Normalize to the mean of each trial (per electrode):
%--------------------------------------------------------------------------
chanmean = mean(squeeze(data(:,:, ntrials)),2);

for nchan = 1:dims(1)
    data(nchan, :, ntrials) = data(nchan, :, ntrials) - chanmean(nchan);

end

end


%%
%--------------------------------------------------------------------------
% Extract band-specific information:
%-------------------------------------------------------------------------- 

%--------------------------------------------------------------------------
% Theta:
%-------------------------------------------------------------------------- 
if band == 1 || strcmpi(band,'theta')

    wa = 4; wb = 7;                     

%--------------------------------------------------------------------------
% Alpha:
%-------------------------------------------------------------------------- 
elseif band == 2 || strcmpi(band,'alpha')

    wa = 8; wb = 12;

%--------------------------------------------------------------------------
% Mu:
%-------------------------------------------------------------------------- 
elseif band == 3 || strcmpi(band,'mu')

    wa = 10; wb = 12;

%--------------------------------------------------------------------------
% Full Beta:
%-------------------------------------------------------------------------- 
elseif band == 4 || strcmpi(band,'full beta')

    wa = 13; wb = 30;

%--------------------------------------------------------------------------
% Low Beta:
%-------------------------------------------------------------------------- 
elseif band == 5 || strcmpi(band,'low beta')

    wa = 13; wb = 18;

%--------------------------------------------------------------------------
% High Beta:
%-------------------------------------------------------------------------- 
elseif band == 6 || strcmpi(band,'high beta')

    wa = 18; wb = 30;
    
%--------------------------------------------------------------------------
% Sensorimotor Beta:
%-------------------------------------------------------------------------- 
elseif band == 7 || strcmpi(band,'sensorimotor beta')

    wa = 18; wb = 22;
    

end

%% 
%--------------------------------------------------------------------------
% Reshape and permute data to optimize filtering:
%-------------------------------------------------------------------------- 
origdims = size(data);

data = reshape(data, [origdims(1), origdims(2) * origdims(3)]);
data = permute(data, [2 1]);

%--------------------------------------------------------------------------
% Build filter:
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Estimate Kaiser parameters for FIR filtering: 
%-------------------------------------------------------------------------- 
[n, Wn, bta, filtype] = kaiserord([wa-2 wa wb wb+2], [0 1 0], [0.01 0.05 0.01], fs);   

%--------------------------------------------------------------------------
% Design FIR filter:
%-------------------------------------------------------------------------- 
b = fir1(n, Wn, filtype, kaiser(n+1, bta),  'noscale');

%--------------------------------------------------------------------------
% Execute filter for desired band:
%-------------------------------------------------------------------------- 
fdata = filtfilt(b, 1, data);

%--------------------------------------------------------------------------
% Undo permutation and then reshaping to regain epochs:
%-------------------------------------------------------------------------- 
data = permute(data, [2 1]);
data = reshape(data, [origdims(1), origdims(2), origdims(3)]);

fdata = permute(fdata, [2 1]);
fdata = reshape(fdata, [origdims(1), origdims(2), origdims(3)]);

%%
%--------------------------------------------------------------------------
% Start loop extract phase representations (per trial/channel):
%-------------------------------------------------------------------------- 
for ntrial = 1:dims(3)
    
    for nchan = 1:dims(1)

%--------------------------------------------------------------------------
% Gabor Representation:
%--------------------------------------------------------------------------
    N = length(fdata(nchan,:,ntrial));
    gabor_rep(nchan,:,ntrial) = hilbert(squeeze(squeeze(fdata(nchan,:,ntrial))),N);

%--------------------------------------------------------------------------
% Phase Construction:
%--------------------------------------------------------------------------
    modulus_rep(nchan,:,ntrial) = abs(squeeze(squeeze(gabor_rep(nchan,:,ntrial))));

%--------------------------------------------------------------------------
% Project to unit circle to obtain phase evolution:
%--------------------------------------------------------------------------
    adata(nchan,:,ntrial) = angle(squeeze(squeeze(gabor_rep(nchan,:,ntrial)))./squeeze(squeeze(modulus_rep(nchan,:,ntrial))));

%--------------------------------------------------------------------------
% Unwrap to get instaneous phase:
%--------------------------------------------------------------------------
    %adata(nchan,:,ntrial) = atan(adata(nchan,:,ntrial));
    %adata(nchan,:,ntrial) = unwrap(adata(nchan,:,ntrial));

%--------------------------------------------------------------------------
% End channel loop:
%--------------------------------------------------------------------------
    end

%--------------------------------------------------------------------------
% End trials loop:
%--------------------------------------------------------------------------
end


%%
pdims = size(adata);

doBase = true;
baseRange = [1 250];

instPLV = zeros(pdims(2), pdims(1), pdims(1));
instiPLV = zeros(pdims(2), pdims(1), pdims(1));
instrPLV = zeros(pdims(2), pdims(1), pdims(1));
bcinstPLV = zeros(pdims(2), pdims(1), pdims(1));

for j = 1:pdims(1)-1
    x = squeeze(adata(j, :, :));
    for k = j+1:pdims(1)
        y = squeeze(adata(k, :, :));
        
            instPLV(:, j, k) = abs(sum(exp(1i*(x - y)), 2))./pdims(3);
            instPLV(:, k, j) = instPLV(:, j, k);

            instiPLV(:, j, k) = abs(imag(sum(exp(1i*(x - y)), 2)))./pdims(3);
            instiPLV(:, k, j) = instiPLV(:, j, k);
            
            instrPLV(:, j, k) = abs(real(sum(exp(1i*(x - y)), 2)))./pdims(3);
            instrPLV(:, k, j) = instrPLV(:, j, k);
            
    end
end

instciPLV =  instiPLV ./ sqrt(1 - (instrPLV.^2));

if true(doBase)
    for j = 1:pdims(1)-1
        for k = j+1:pdims(1)
            mubl = mean(instPLV(baseRange(1):baseRange(2), j, k),1);
            sdbl = std(instPLV(baseRange(1):baseRange(2), j, k),1);
            bcinstPLV(:, j, k) = (instPLV(:, j, k) - mubl)/sdbl;
            bcinstPLV(:, k, j) = bcinstPLV(:, j, k);

            mubl = mean(instiPLV(baseRange(1):baseRange(2), j, k),1);
            sdbl = std(instiPLV(baseRange(1):baseRange(2), j, k),1);
            bcinstiPLV(:, j, k) = (instiPLV(:, j, k) - mubl)/sdbl;
            bcinstiPLV(:, k, j) = bcinstiPLV(:, j, k);
            
            mubl = mean(instciPLV(baseRange(1):baseRange(2), j, k),1);
            sdbl = std(instciPLV(baseRange(1):baseRange(2), j, k),1);
            bcinstciPLV(:, j, k) = (instciPLV(:, j, k) - mubl)/sdbl;
            bcinstciPLV(:, k, j) = bcinstciPLV(:, j, k);
            
        end
    end

end
            
            
instPLV = squeeze(instPLV);



if true(doBase)
    bcinstPLV = squeeze(bcinstPLV);
    
end


end