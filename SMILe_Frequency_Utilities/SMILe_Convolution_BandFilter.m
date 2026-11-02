function [fdata, as] = SMILe_Convolution_BandFilter(data, fs, lowerBound, upperBound, transitionWidth)
%--------------------------------------------------------------------------
% Filters a data set using convolution filtering. The user defines the 
% lower and upper bounds of the frequency range and the "sharpness" of the 
% filter. 
%
% Input:
% data: Data to be filtered. Typically a 3-D matrix of chan x time x trials
%       Can also be a 2-D matrix or 1-D vector.
%
% fs:   The sampling rate (in Hz) at which the input data were collected.
%
% lowerBound: Lower boundary of the frequency range of interest
%
% upperBound: Upper boundary of the frequency range of interest
% 
% transitionWidth: Sets the "sharpness" of the filter.
%
% Output:
% fdata: The filtered data, averaged across the individual frequencies
%       within the band of interest.
%
% as: The analytic signal averaged across the individual frequencies
%       within the band of interest.
%
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
        chanmean = mean(squeeze(data(:,:, ntrials)),2);

            for nchan = 1:dims(1)
                data(nchan, :, ntrials) = data(nchan, :, ntrials) - chanmean(nchan);

            end

        end
       
%--------------------------------------------------------------------------
% Reshape and permute data to optimize filtering:
%--------------------------------------------------------------------------
        origdims = size(data);

        data = reshape(data, [origdims(1), origdims(2) * origdims(3)]);
               
        if dims(1) > dims(2)
            data = permute(data,[2 1]);
           
        end
       
        %data = permute(data, [2 1]);
       
    case 2
       
        if dims(1) > dims(2)
            data = permute(data,[2 1]);
           
        end
       
        chanmean = mean(squeeze(data(:,:)), 2);
       
            for nchan = 1:dims(1)
                data(nchan, :) = data(nchan, :) - chanmean(nchan);

            end
           
        origdims = size(data);
        %data = permute(data, [2 1]);
           
    case 1
           
        chanmean = mean(squeeze(data), 2);
        data = data - chanmean;
           
        origdims = size(data);
           
            if dims(1) > dims(2)
            data = permute(data,[2 1]);
           
            end
           
            %data = permute(data, [2 1]);
           
end
       

%%
%--------------------------------------------------------------------------
% Build filter:
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Specify Nyquist freuqency
%--------------------------------------------------------------------------
nyquist = fs/2;

%--------------------------------------------------------------------------
% Design filter parameters that approximately meets the specifications:
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% A typical filter order: 3 times the ratio of the sampling rate to the 
% lower edge of the frequency band.
%--------------------------------------------------------------------------
filt_order = round(3*(fs/lowerBound)); 

%--------------------------------------------------------------------------
% Set frequency vector as fraction of Nyquist:
%--------------------------------------------------------------------------
ffrequencies  = [ 0 (1-transitionWidth)*lowerBound lowerBound upperBound (1+transitionWidth)*upperBound nyquist ]/nyquist;

%--------------------------------------------------------------------------
% Shape of filter must be the same number of elements as frequency vector:
%--------------------------------------------------------------------------
idealresponse = [ 0 0 1 1 0 0 ];

%--------------------------------------------------------------------------
% Compute filter weights:
%--------------------------------------------------------------------------
filterweights = firls(filt_order,ffrequencies,idealresponse);
 
%--------------------------------------------------------------------------
% Prep fft and ifft settings:
%--------------------------------------------------------------------------
fftDims = size(data,2);

%--------------------------------------------------------------------------
% Get length of filter kernel:
%--------------------------------------------------------------------------
nKernel = length(filterweights);

%--------------------------------------------------------------------------
% Need to set up to cut output of the convolution:
%--------------------------------------------------------------------------
halfKernel = floor(nKernel/2);

%--------------------------------------------------------------------------
% Optomize speed of calculations:
%--------------------------------------------------------------------------
nConv = fftDims + nKernel - 1;
nConvPow2 = pow2(nextpow2(nConv));

%--------------------------------------------------------------------------
% Calculate Fourier of the filter:
%--------------------------------------------------------------------------
fftKernel = fft(filterweights, nConvPow2);
fftKernel = fftKernel / max(fftKernel);

%--------------------------------------------------------------------------
% Allocate memory for output:
%--------------------------------------------------------------------------
fdata = zeros(size(data));
as = zeros(size(data));

%--------------------------------------------------------------------------
% Start progress bar:
%--------------------------------------------------------------------------
h = waitbar(0,'Starting...');

%--------------------------------------------------------------------------
% Use fft() and ifft() (manual convolution) and the filter weights to 
% filter the signal:
%--------------------------------------------------------------------------
for nchan = 1:origdims(1)
        waitbar(nchan/origdims(1), h, sprintf('Filtering channel: %6.0f ', nchan))

%--------------------------------------------------------------------------
% Fourier transform of the signal:
%--------------------------------------------------------------------------        
        fftData = fft(data(nchan,:), nConvPow2);

%--------------------------------------------------------------------------
% Convolve the filter weights and data to obtain the analytic signal:
%--------------------------------------------------------------------------        
        convResult = (fftKernel .* fftData);
        convResult = convResult(1:nConv);
        as(nchan, :) = convResult(halfKernel+1 : end-halfKernel);

%--------------------------------------------------------------------------
% Return filtered signal to the time domain:
%--------------------------------------------------------------------------                
        filtResult = ifft(convResult);
        filtResult = filtResult(1:nConv);
        filtResult = filtResult(halfKernel+1 : end-halfKernel);
        fdata(nchan, :) = real(filtResult)*2;

%--------------------------------------------------------------------------
% Alternatively, use the built-in convolution function:
%--------------------------------------------------------------------------
        %fdata(nchan, :) = conv(data(nchan,:),filterweights,'same');
    
end

%--------------------------------------------------------------------------
% Kill progress bar:
%--------------------------------------------------------------------------
close(h)

%--------------------------------------------------------------------------
% Undo reshape to regain epochs if those came in:
%--------------------------------------------------------------------------
if numdims == 3
    fdata = reshape(fdata, [origdims(1), origdims(2), origdims(3)]);
    as = reshape(as, [origdims(1), origdims(2), origdims(3)]);    

end

%--------------------------------------------------------------------------
% For visualization
%--------------------------------------------------------------------------
% figure;
% hold on
% plot(squeeze(fdata(17,:,5)))
% plot(abs(squeeze(fdata(17,:,5))),'r')
% %[yupper ylower] = envelope(abs(squeeze(fdata(17,:,5))),20,'peak');
% [env] = SMILe_Freq_Envelope(abs(squeeze(fdata(17,:,5))),10,'peak');
% plot(env.upper,'k')

%--------------------------------------------------------------------------
% End of function
%--------------------------------------------------------------------------
end
