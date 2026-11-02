function [surr] = SMILe_Surrogates_RandPhase_RandEpochs(datain, nreps)
%--------------------------------------------------------------------------
%
% This script generates phase randomized surrogate data that preserves
% spatial and temporal correlation of the original data.
%
% Based on Prichard D, Theiler J. Generating surrogate data for time series
% with several simultaneously measured variables. Physical review letters.
% 1994 Aug 15;73(7):951.
%
% Input:
% datain:   An EEGLAB-type data set arranged in channels x time x epochs
%
% nreps:    Number of complete repetitions of surrogate eopchs to be
%           created. If the original data have multiple epochs, the data 
%           will be averaged across epochs, nreps surrogate epochs will be 
%           created.
%
% Output:
% surr:     channels x time x N surrogate epochs
%
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Check the original dimensions, and reshape for convenience. Mark the data
% for re-permuting if altered:
%--------------------------------------------------------------------------
origdims = size(datain);

if numel(origdims) == 3
    datain = squeeze(mean(datain,3));

end

if origdims(2) > origdims(1)
    datain = permute(datain, [2 1]);
    repermute = true;
   
end

%--------------------------------------------------------------------------
% New dimensions:
%--------------------------------------------------------------------------
dims = size(datain);

%--------------------------------------------------------------------------
% Enforce an even time range:
%--------------------------------------------------------------------------
etime = round(dims(1)/2)*2;

%--------------------------------------------------------------------------
% Start the randomization for surrogates with a null data set:
%--------------------------------------------------------------------------
surr = [];

%--------------------------------------------------------------------------
% The fft call will add a zero at the end if the length is uneven:
%--------------------------------------------------------------------------
Xfft = fft(datain, etime);
       
%--------------------------------------------------------------------------
% Get the original amplitude:
%--------------------------------------------------------------------------
amp = abs(Xfft(1:etime/2+1, :));
       
%--------------------------------------------------------------------------
% Get the original phase:
%--------------------------------------------------------------------------
phi = angle(Xfft(1:etime/2+1, :));

%--------------------------------------------------------------------------
% Start loop around the number of surrogate data sets:
%--------------------------------------------------------------------------
for counter = 1:nreps
       
%--------------------------------------------------------------------------
% Calculate the random phase to add:
%--------------------------------------------------------------------------
    phir = 4*acos(0)*rand(etime/2-1,1)-2*acos(0);
       
%--------------------------------------------------------------------------
% Apply random phases to the original data (complex plane):
%--------------------------------------------------------------------------
    tmp(2:etime/2,:) = amp(2:etime/2, :).*exp(sqrt(-1)*(phi(2:etime/2, :) + repmat(phir, 1, dims(2))));
       
%--------------------------------------------------------------------------
% Take the inverse FFT to bring evertying back to the real plane:
%--------------------------------------------------------------------------
    tmp = ifft([Xfft(1, :); tmp(2:etime/2, :); Xfft(etime/2+1, :); conj(tmp(etime/2:-1:2, :))]);
       
%--------------------------------------------------------------------------
% Collect epochs:
%--------------------------------------------------------------------------
    randdat = tmp(1:dims(1), :);
       
%--------------------------------------------------------------------------
% Concatenate data sets across the number of repetitions:
%--------------------------------------------------------------------------
    surr = cat(3, surr, randdat);

%--------------------------------------------------------------------------
% End loop for surrogate repetitions:
%--------------------------------------------------------------------------
end

%--------------------------------------------------------------------------
% Mean normalize all of the surrogates:
%--------------------------------------------------------------------------
bl = mean(surr(:,:), 2);
surr = bsxfun(@minus, bl, surr);

%--------------------------------------------------------------------------
% Reshape to original format, if needed:
%--------------------------------------------------------------------------
if true(repermute)
    surr = permute(surr, [2 1 3]);
   
end

%--------------------------------------------------------------------------
% End of function
%--------------------------------------------------------------------------
end