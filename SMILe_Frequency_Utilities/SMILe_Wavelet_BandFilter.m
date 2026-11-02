function [fdata, pow, as, cs] = SMILe_Wavelet_BandFilter(data,Fs,fWavelet, doCS)
%--------------------------------------------------------------------------
% Filters a data set using wavelet-based filtering. The user defines the 
% wavelet center frequencies. 
%
% Input:
% data: Data to be filtered. Typically a 3-D matrix of chan x time x trials
%       Can also be a 2-D matrix or 1-D vector.
%
% Fs:   The sampling rate (in Hz) at which the input data were collected.
% 
% fWavelet: Wavelet center frequencies. In the frequency domain, 
%       the wavelets consist of a Gaussian function of frequency with 
%       center frequency represented by the peak of the Gaussian. 
%       Center frequencies (fWavelet) should consist of a single frequency
%       or multiple frequencies (e.g. [10], or [8:12]). Frequencies 
%       are in Hz.
%
% Output:
% fdata: The filtered data, averaged across the individual frequencies
%       within the band of interest.
%
% pow: Time-resolved absolute power (non-normalized) averaged across the
%       individual frequencies within the band of interest.
%
% as: The analytic signal averaged across the individual frequencies
%       within the band of interest.
%
% cs: The cross spectrum taken from the analytic signal, averaged across 
%       the individual frequencies within the band of interest.
%%
%--------------------------------------------------------------------------
% Confirm double precision:
%--------------------------------------------------------------------------
data = double(data);

%--------------------------------------------------------------------------
% Get dimensions of the data, use switch/case for dimensions:
%--------------------------------------------------------------------------
dims = size(data);
numdims = numel(dims);
if min(dims) == 1
    numdims = 1;
   
end

switch numdims
   
    case 3
%--------------------------------------------------------------------------
% Loop across trials:
%--------------------------------------------------------------------------
        for ntrials = 1:dims(3)
   
%--------------------------------------------------------------------------
% Normalize to the mean of each trial (per electrode or for a single
% signal):
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Also reshape and permute data to optimize filtering:
%--------------------------------------------------------------------------
        chanmean = mean(squeeze(data(:,:, ntrials)),2);

            for nchan = 1:dims(1)
                data(nchan, :, ntrials) = data(nchan, :, ntrials) - chanmean(nchan);

            end

        end
       
        origdims = size(data);

        data = reshape(data, [origdims(1), origdims(2) * origdims(3)]);
               
        if dims(1) > dims(2)
            data = permute(data,[2 1]);
           
        end
       
    case 2
       
        if dims(1) > dims(2)
            data = permute(data,[2 1]);
           
        end
       
        chanmean = mean(squeeze(data(:,:)), 2);
       
            for nchan = 1:dims(1)
                data(nchan, :) = data(nchan, :) - chanmean(nchan);

            end
           
        origdims = size(data);
           
    case 1
           
        chanmean = mean(squeeze(data), 2);
        data = data - chanmean;
           
        origdims = size(data);
           
            if dims(1) > dims(2)
            data = permute(data,[2 1]);
           
            end
           
end

%--------------------------------------------------------------------------
% Begin wavelet creation:
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Get frequency characteristics:
%--------------------------------------------------------------------------
minFreq = fWavelet(1);
maxFreq = fWavelet(end);
numFreq = numel(fWavelet);
frex = linspace(minFreq,maxFreq,numFreq);

%--------------------------------------------------------------------------
% Parameterize the range of cycles:
%--------------------------------------------------------------------------
rangeCycles = [ 4 10 ];

%--------------------------------------------------------------------------
% Logarithmic or linear spacing?
%--------------------------------------------------------------------------
%wavecycles = logspace(log10(rangeCycles(1)),log10(rangeCycles(end)),numFreq) ./ (2*pi*frex);
wavecycles = linspace(rangeCycles(1), rangeCycles(end), numFreq) ./ (2*pi*frex);

%--------------------------------------------------------------------------
% Wavelet duration/length:
%--------------------------------------------------------------------------
wavtime = -2:1/Fs:2;
halfWave = (length(wavtime)-1)/2;

%--------------------------------------------------------------------------
% Set FFT parameters:
%--------------------------------------------------------------------------
nWave = length(wavtime);

if numdims == 3
    nData = dims(2) * dims(3); 

else
    nData = dims(2);

end

nConv = nWave + nData - 1;

%--------------------------------------------------------------------------
% Start waitbar:
%--------------------------------------------------------------------------
multiWaitbar( 'Frequencies', 'Color', [0.2 0.9 0.3] );
multiWaitbar( 'Channels', 'Color', [0.8 0.0 0.1] );
if doCS == 1
    multiWaitbar( 'Cross Spectrum', 'Color', [0.1 0.0 0.8] );
end

%--------------------------------------------------------------------------
% Allocate memory:
%--------------------------------------------------------------------------
if numdims == 3
    f = zeros(dims(1),dims(2)*dims(3),numel(fWavelet));
    xas = zeros(dims(1),dims(2)*dims(3),numel(fWavelet));
    tfpow = zeros(dims(1),dims(2)*dims(3),numel(fWavelet));

else
    f = zeros(dims(1),dims(2),numel(fWavelet));
    xas = zeros(dims(1),dims(2),numel(fWavelet));
    tfpow = zeros(dims(1),dims(2),numel(fWavelet));

end

%--------------------------------------------------------------------------
% Loop over frequencies:
%--------------------------------------------------------------------------
for numfreq = 1:length(frex)
    
%--------------------------------------------------------------------------
% Update waitbar:
%--------------------------------------------------------------------------
    multiWaitbar( 'Frequencies', numfreq/length(frex));

%--------------------------------------------------------------------------
% Create wavelet and get its FFT:
%--------------------------------------------------------------------------
    wavelet  = exp(2*1i*pi*frex(numfreq).*wavtime) .* exp(-wavtime.^2./(2*wavecycles(numfreq)^2));
    waveletX = fft(wavelet,nConv);
    waveletX = waveletX ./ max(waveletX);
    
    for nchan = 1:size(data,1)

%--------------------------------------------------------------------------
% Update waitbar:
%--------------------------------------------------------------------------
    multiWaitbar( 'Channels', nchan/dims(1) );

%--------------------------------------------------------------------------
% Select channel data and get the Fourier transform:
%--------------------------------------------------------------------------
        dataX = fft(data(nchan,:),nConv);

%--------------------------------------------------------------------------
% Run the convolution:
%--------------------------------------------------------------------------
        convResult = ifft(waveletX .* dataX);
        convResult = convResult(halfWave+1:end-halfWave);
        f(nchan,:,numfreq) = 2.*real(convResult);
        xas(nchan,:,numfreq) = convResult;
        tfpow(nchan,:,numfreq) = abs(convResult).^2;

    end

end

%--------------------------------------------------------------------------
% Collect output:
%--------------------------------------------------------------------------
fdata = mean(f,3);
as = mean(xas,3);
pow = mean(tfpow,3);

%--------------------------------------------------------------------------
% Compute the cross spectrum from the analytic signal:
%--------------------------------------------------------------------------
if exist('doCS', 'var') && ~isempty(doCS) && doCS > 0
    cs = zeros(nchan, nchan, size(as,2));
    
    for chani = 1:nchan
        multiWaitbar( 'Cross Spectrum', chani/nchan);
        
        for chanj = chani:nchan
            cs(chani, chanj, :,:) = squeeze(as(chani,:,:) .* conj(as(chanj,:,:)));
            
        end
        
    end
    
    cs = squeeze(cs);

elseif ~exist('doCS', 'var') && isempty(doCS) || doCS == 0
    cs = [];
    
end

%--------------------------------------------------------------------------
% Regain epochs if those came in:
%--------------------------------------------------------------------------
if numdims == 3
    fdata = reshape(fdata, [origdims(1), origdims(2), origdims(3)]);
    as = reshape(as, [origdims(1), origdims(2), origdims(3)]);
    pow = reshape(pow, [origdims(1), origdims(2), origdims(3)]);
    if doCS == 1
        cs = reshape(cs, [origdims(1), origdims(1), origdims(2), origdims(3)]);
    end
    
end

%--------------------------------------------------------------------------
% Kill the waitbars:
%--------------------------------------------------------------------------
multiWaitbar('CloseAll');

    
end