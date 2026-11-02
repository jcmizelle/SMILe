function [fdata] = SMILe_FIR_BandFilter(data, fs, lowerBound, upperBound, transitionWidth)
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
       

%%

%--------------------------------------------------------------------------
% Build filter:
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Estimate parameters for FIR filtering. This estimates the minimum
% FIR filter order that will approximately meet the given set of
% specifications:
%--------------------------------------------------------------------------
nyquist = fs/2;

%--------------------------------------------------------------------------
% A typical filter order: 3 times the ratio of the sampling rate to the 
% lower edge of the frequency band.
%--------------------------------------------------------------------------
filterOrder = round(3*(fs/lowerBound));

%--------------------------------------------------------------------------
% Set frequency vector as fraction of Nyquist:
%--------------------------------------------------------------------------
frex  = [ 0 (1-transitionWidth)*lowerBound lowerBound upperBound (1+transitionWidth)*upperBound nyquist ]/nyquist;

%--------------------------------------------------------------------------
% Shape of filter must be the same number of elements as frequency vector:
%--------------------------------------------------------------------------
idealResponse = [ 0 0 1 1 0 0 ];

%--------------------------------------------------------------------------
% Compute filter weights:
%--------------------------------------------------------------------------
filterweights = firls(filterOrder, frex, idealResponse);

%--------------------------------------------------------------------------
% Performs zero-phase digital filtering by processing the input data in
% both the forward and reverse directions.
%
% filtfilt minimizes start-up and ending transients by matching initial
% conditions, and works for both real and complex inputs.
%
% Execute filter for desired band, zeropad if necessary:
%--------------------------------------------------------------------------
fdata = zeros(size(data));
for nchan = 1:origdims(1)
    fdata(nchan, :) = filtfilt(filterweights,1,data(nchan,:));
   
end

%--------------------------------------------------------------------------
% Undo permutation and then reshape to regain epochs if those came in:
%--------------------------------------------------------------------------

if numdims == 3
    fdata = reshape(fdata, [origdims(1), origdims(2), origdims(3)]);

end

%--------------------------------------------------------------------------
% End for function:
%--------------------------------------------------------------------------


%% For visualization
% figure;
% hold on
% plot(squeeze(fdata(17,:,5)))
% plot(abs(squeeze(fdata(17,:,5))),'r')
% %[yupper ylower] = envelope(abs(squeeze(fdata(17,:,5))),20,'peak');
% [env] = SMILe_Freq_Envelope(abs(squeeze(fdata(17,:,5))),10,'peak');
% plot(env.upper,'m')

end