function [surr] = SMILe_Surrogates_RandPhase(datain,nsurr)
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
%%
%--------------------------------------------------------------------------
% Check the original dimensions, and reshape for convenience. Mark the data
% for re-permuting if altered:
%--------------------------------------------------------------------------
origdims = size(datain);
numchan = origdims(1);
numdims = numel(origdims);


if numdims == 3
    datain = squeeze(mean(datain,3));

end

%--------------------------------------------------------------------------
% New dimensions:
%--------------------------------------------------------------------------
dims = size(datain);

%--------------------------------------------------------------------------
% Allocate memory for temporary storage of surrogates:
%--------------------------------------------------------------------------
outsurr = zeros([length(datain) nsurr]);

%--------------------------------------------------------------------------
% Allocate memory for output:
%--------------------------------------------------------------------------
surr = zeros([dims nsurr]);

%--------------------------------------------------------------------------
% Start progress bar:
%--------------------------------------------------------------------------
%wh = waitbar(0,'Starting...');

%multiWaitbar( 'Channels', 'Color', [0.8 0.0 0.1] );
%multiWaitbar( 'Surrogates', 'Color', [0.2 0.9 0.3] );

%--------------------------------------------------------------------------
% Start loop across channels:
%--------------------------------------------------------------------------
for nchan = 1:numchan
    %waitbar(nchan/numchan, wh, sprintf('Creating surrogates for channel: %6.0f ', nchan))
%--------------------------------------------------------------------------
% Update waitbar:
%--------------------------------------------------------------------------
    %multiWaitbar( 'Channels', nchan/numchan );

%--------------------------------------------------------------------------
% Pull data channel-by-channel:
%--------------------------------------------------------------------------
    sig = datain(nchan,:);
    sig = transpose(sig);
    sigleng = length(sig);

%--------------------------------------------------------------------------
% Start loop around the number of surrogate data sets:
%--------------------------------------------------------------------------
    %for counter = 1:nsurr
%--------------------------------------------------------------------------
% Fourier transform of the signal:
%--------------------------------------------------------------------------
        y = fft(sig);

%--------------------------------------------------------------------------
% Magnitude of the spectrum:
%--------------------------------------------------------------------------
        m = abs(y);
%--------------------------------------------------------------------------
% Calculate phase angles:
%--------------------------------------------------------------------------
        p = angle(y);
%--------------------------------------------------------------------------
% Calculate the imaginary unit:
%--------------------------------------------------------------------------
        i = sqrt(-1);

%--------------------------------------------------------------------------
% Work on half of the data points:
%--------------------------------------------------------------------------
        h = floor(sigleng/2);

%--------------------------------------------------------------------------
% Loop across surrogates:
%--------------------------------------------------------------------------
        for  j = 1:nsurr
            %multiWaitbar( 'Surrogates', j/nsurr );
%--------------------------------------------------------------------------
% Randomize phases:
%--------------------------------------------------------------------------            
           if rem(sigleng,2) == 0
              p1 = rand(h-1,1)*2*pi;
              p(2:sigleng) = [p1' p(h+1) -flipud(p1)'];
%--------------------------------------------------------------------------
% Adjust magnitudes:
%--------------------------------------------------------------------------              
              m = [m(1:h+1);flipud(m(2:h))];
              
           else
              p1 = rand(h,1)*2*pi;
              p(2:sigleng) = [p1 -flipud(p1)];
              
           end
%--------------------------------------------------------------------------
% Back to complex numbers:
%--------------------------------------------------------------------------
           outsurr(:,j) = m.*exp(i*p);
%--------------------------------------------------------------------------
% Return the time series (phase randomized surrogates):
%--------------------------------------------------------------------------           
           outsurr(:,j) = real(ifft(outsurr(:,j)));
           
        end


%--------------------------------------------------------------------------
% Collect output:
%--------------------------------------------------------------------------
surr(nchan,:,:) = outsurr;

%--------------------------------------------------------------------------
% Refresh chansurrs:
%--------------------------------------------------------------------------
outsurr = zeros([length(datain) nsurr]);

end

%--------------------------------------------------------------------------
% Kill progress bar:
%--------------------------------------------------------------------------
%close(wh)
%multiWaitbar('CloseAll');

%--------------------------------------------------------------------------
% End for function:
%--------------------------------------------------------------------------
end