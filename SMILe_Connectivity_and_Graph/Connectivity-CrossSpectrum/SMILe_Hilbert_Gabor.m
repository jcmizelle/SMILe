function [gabor_as, gabor_mod, phaseAngle, instantPhase, sigpow] = SMILe_Hilbert_Gabor(fdata)
%%
%--------------------------------------------------------------------------
%
% The Hilbert transform turns a real-valued signal into a complex-valued 
% signal, also known as the analytic signal, a complex signal defined as: 
% x_a(t) = x(t) + iHT{x(t)} = Re(x_a(t)) + Im(x_a(t))
% which means turning a real-valued 1-D time series into a complex-valued 
% 2-D time series with both real and imaginary components. This is also 
% known as the Gabor representation of the signal.
%
% The ultimate outcome is similar to that coming from the Morelet Wavelet
% decomposition that yields time and frequency information. The wavelet
% transform is superior for broad-spectrum data - such as EEG data that
% have not been band-pass filtered for a specific frequency.
% 
% Once the complex analytic signal has been obtained, there are other
% measures that can be derived, such as the instantaneous phase angle (for
% connectivity analysis) and instantaneous power.
%
% The instantaneous power is not computed on the filtered signal, but the 
% analytic signal. Namely, it is the squared magnitude of the complex 
% vector at an instant in time. Similarly, the instantaneous phase is the 
% phase angle between the vector and the real axis (0 deg) in the complex 
% plane. Mathematically, this refers to the complex sinusoid, i.e., 
% x(t) = Ae^(iwt) in Euler form, or x(t) = A(cos(wt)+isin(wt)) in 
% rectangular form.
%
%--------------------------------------------------------------------------
%
% Input:
% - fdata: Must be band-pass filtered for a frequency of interest. See 
%   SMILe_BandFilter for how to obtain this.
%
% Output:
% - gabor_as: complex valued analytic signal
%
% - phaseAngle: instantaneous phase angle
%
% - instantPhase: the instantaneous phase of the signal
%
% - sigpow: the instantaneous power of the signal
%
%--------------------------------------------------------------------------

%%
%--------------------------------------------------------------------------
% Grab information related to the size of incoming data set:
%--------------------------------------------------------------------------    
dims = size(fdata);
numdims = numel(dims);

switch numdims
    
%--------------------------------------------------------------------------
% chans x time x trials (or participants):
%--------------------------------------------------------------------------            
    case 3
%--------------------------------------------------------------------------
% Preallocate memory:
%--------------------------------------------------------------------------                    
        gabor_as = zeros(dims);
        gabor_mod = zeros(dims);
        phaseAngle = zeros(dims);
        instantPhase = zeros(dims);
        sigpow = zeros(dims);

%--------------------------------------------------------------------------
% Reshape the data to improve "filtering" performance at the edges of the 
% analytic signal:
%--------------------------------------------------------------------------          
        fdata = reshape(fdata, [dims(1), dims(2) * dims(3)]);
        
%--------------------------------------------------------------------------
% Start loop for trials:
%--------------------------------------------------------------------------                  
        for nchan = 1:dims(1)
            data = squeeze(fdata(nchan,:));
            npts = size(data,2);
        
%--------------------------------------------------------------------------
% Calculate the Fast Fourier Transform of the (real) signal:
%--------------------------------------------------------------------------        
            fqdomain = fft(data);
        
%--------------------------------------------------------------------------
% Create a copy that is multiplied by the complex operator:
%--------------------------------------------------------------------------
            complexfq = 1i*fqdomain;

%--------------------------------------------------------------------------
% Find the indices of positive and negative frequencies:
%--------------------------------------------------------------------------
            posFq = 2:floor(npts/2) + mod(npts,2);
            negFq = ceil(npts/2)+1 + ~mod(npts,2):npts;

%--------------------------------------------------------------------------
% Rotate the Fourier Coefficients:
% 
% This works by computing the iAsin(2pft) component (i.e., the phase 
% quadrature).
%
% Positive frequencies are rotated counter-clockwise; negative frequencies 
% are rotated clockwise).
%
% Can substitute this as an alternative and slightly faster method:
% f(posF) = f(posF)*2;
% f(negF) = f(negF)*0;
%--------------------------------------------------------------------------
            fqdomain(posFq) = fqdomain(posFq) + -1i*complexfq(posFq);
            fqdomain(negFq) = fqdomain(negFq) +  1i*complexfq(negFq);

%--------------------------------------------------------------------------
% Take inverse FFT to compute the Hilbert transform, which provides the
% analytic signal to use for subsequent analyses. This is also known as the
% Gabor representaton
%--------------------------------------------------------------------------
            gabor_as(nchan,:) = ifft(fqdomain);
        
%--------------------------------------------------------------------------
% Calculate the modulus of the Gabor representation. This also represents 
% the amplitude envelope:
%--------------------------------------------------------------------------
            gabor_mod(nchan,:) = abs(gabor_as(nchan,:));
        
%--------------------------------------------------------------------------
% Phase Construction:
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Project to unit circle to obtain phase evolution:
%--------------------------------------------------------------------------
            %phaseAngle = angle(gabor_as);
            phaseAngle(nchan,:) = angle(gabor_as(nchan,:)./gabor_mod(nchan,:));

%--------------------------------------------------------------------------
% Unwrap to get instaneous phase:
%--------------------------------------------------------------------------
            instantPhase(nchan,:) = unwrap(phaseAngle(nchan,:));

%--------------------------------------------------------------------------
% Instantaneous power (for the specific frequency band brought into the 
% function) can also be obtained:
%--------------------------------------------------------------------------
            sigpow(nchan,:) = (abs(gabor_as(nchan,:))).^2;
            sigpow(nchan,:) = 10*log10(sigpow(nchan,:));

%--------------------------------------------------------------------------
% End channel loop:
%--------------------------------------------------------------------------
        end
    
%--------------------------------------------------------------------------
% Recapture original shape of the data:
%--------------------------------------------------------------------------
    gabor_as = reshape(gabor_as, [dims(1), dims(2), dims(3)]);
    gabor_mod = reshape(gabor_mod, [dims(1), dims(2), dims(3)]);
    phaseAngle = reshape(phaseAngle, [dims(1), dims(2), dims(3)]);
    instantPhase = reshape(instantPhase, [dims(1), dims(2), dims(3)]);
    sigpow = reshape(sigpow, [dims(1), dims(2), dims(3)]);

%--------------------------------------------------------------------------
% chans x time:
%--------------------------------------------------------------------------        
    case 2
%--------------------------------------------------------------------------
% Preallocate memory:
%--------------------------------------------------------------------------                    
        gabor_as = zeros(dims);
        gabor_mod = zeros(dims);
        phaseAngle = zeros(dims);
        instantPhase = zeros(dims);
        sigpow = zeros(dims);
        
        for nchan = 1:dims(1)
            data = squeeze(fdata(nchan,:));
            npts = size(data,2);
        
%--------------------------------------------------------------------------
% Calculate the Fast Fourier Transform of the (real) signal:
%--------------------------------------------------------------------------        
            fqdomain = fft(data);
        
%--------------------------------------------------------------------------
% Create a copy that is multiplied by the complex operator:
%--------------------------------------------------------------------------
            complexfq = 1i*fqdomain;

%--------------------------------------------------------------------------
% Find the indices of positive and negative frequencies:
%--------------------------------------------------------------------------
            posFq = 2:floor(npts/2) + mod(npts,2);
            negFq = ceil(npts/2)+1 + ~mod(npts,2):npts;

%--------------------------------------------------------------------------
% Rotate the Fourier Coefficients:
% 
% This works by computing the iAsin(2pft) component (i.e., the phase 
% quadrature).
%
% Positive frequencies are rotated counter-clockwise; negative frequencies 
% are rotated clockwise).
%
% Can substitute this as an alternative and slightly faster method:
% f(posF) = f(posF)*2;
% f(negF) = f(negF)*0;
%--------------------------------------------------------------------------
            fqdomain(posFq) = fqdomain(posFq) + -1i*complexfq(posFq);
            fqdomain(negFq) = fqdomain(negFq) +  1i*complexfq(negFq);

%--------------------------------------------------------------------------
% Take inverse FFT to compute the Hilbert transform, which provides the
% analytic signal to use for subsequent analyses. This is also known as the
% Gabor representaton
%--------------------------------------------------------------------------
            gabor_as(nchan,:) = ifft(fqdomain);
        
%--------------------------------------------------------------------------
% Calculate the modulus of the Gabor representation. This also represents 
% the amplitude envelope:
%--------------------------------------------------------------------------
            gabor_mod(nchan,:) = abs(gabor_as(nchan,:));

%--------------------------------------------------------------------------
% Project to unit circle to obtain phase evolution:
%--------------------------------------------------------------------------
            %phaseAngle = angle(gabor_as);
            phaseAngle(nchan,:) = angle(gabor_as(nchan,:)./gabor_mod(nchan,:));

%--------------------------------------------------------------------------
% Unwrap to get instaneous phase:
%--------------------------------------------------------------------------
            instantPhase(nchan,:) = unwrap(phaseAngle(nchan,:));

%--------------------------------------------------------------------------
% Instantaneous power (for the specific frequency band brought into the 
% function) can also be obtained:
%--------------------------------------------------------------------------
            sigpow(nchan,:) = (abs(gabor_as(nchan,:))).^2;

%--------------------------------------------------------------------------
% End channel loop:
%--------------------------------------------------------------------------
        end

%--------------------------------------------------------------------------
% single vector:
%--------------------------------------------------------------------------        
    case 1

                data = squeeze(fdata);
                npts = size(data,2);
        
%--------------------------------------------------------------------------
% Calculate the Fast Fourier Transform of the (real) signal:
%--------------------------------------------------------------------------        
                fqdomain = fft(data);
        
%--------------------------------------------------------------------------
% Create a copy that is multiplied by the complex operator:
%--------------------------------------------------------------------------
                complexfq = 1i*fqdomain;

%--------------------------------------------------------------------------
% Find the indices of positive and negative frequencies:
%--------------------------------------------------------------------------
                posFq = 2:floor(npts/2) + mod(npts,2);
                negFq = ceil(npts/2)+1 + ~mod(npts,2):npts;

%--------------------------------------------------------------------------
% Rotate the Fourier Coefficients:
% 
% This works by computing the iAsin(2pft) component (i.e., the phase 
% quadrature).
%
% Positive frequencies are rotated counter-clockwise; negative frequencies 
% are rotated clockwise).
%
% Can substitute this as an alternative and slightly faster method:
% f(posF) = f(posF)*2;
% f(negF) = f(negF)*0;
%--------------------------------------------------------------------------
                fqdomain(posFq) = fqdomain(posFq) + -1i*complexfq(posFq);
                fqdomain(negFq) = fqdomain(negFq) +  1i*complexfq(negFq);

%--------------------------------------------------------------------------
% Take inverse FFT to compute the Hilbert transform, which provides the
% analytic signal to use for subsequent analyses. This is also known as the
% Gabor representaton
%--------------------------------------------------------------------------
                gabor_as = ifft(fqdomain);
        
%--------------------------------------------------------------------------
% Calculate the modulus of the Gabor representation. This also represents 
% the amplitude envelope:
%--------------------------------------------------------------------------
                gabor_mod = abs(gabor_as);

%--------------------------------------------------------------------------
% Project to unit circle to obtain phase evolution:
%--------------------------------------------------------------------------
                %phaseAngle = angle(gabor_as);
                phaseAngle = angle(gabor_as./gabor_mod);

%--------------------------------------------------------------------------
% Unwrap to get instaneous phase:
%--------------------------------------------------------------------------
                instantPhase = unwrap(phaseAngle);

%--------------------------------------------------------------------------
% Instantaneous power (for the specific frequency band brought into the 
% function) can also be obtained:
%--------------------------------------------------------------------------
                sigpow = (abs(gabor_as)).^2;

%--------------------------------------------------------------------------
% End for switch/case: 
%--------------------------------------------------------------------------
end

%--------------------------------------------------------------------------
% End for function:
%--------------------------------------------------------------------------
end

% figure(3), clf
% plot3(EEG.times,real(hilbertx),imag(hilbertx))
% xlabel('Time (ms)'), ylabel('Real part'), zlabel('Imaginary part')
% axis tight
% rotate3d


