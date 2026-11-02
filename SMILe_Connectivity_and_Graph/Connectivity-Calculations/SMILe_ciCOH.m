function [icoh, cicoh] = SMILe_ciCOH(cs)
%--------------------------------------------------------------------------
%
% Calculates imaginary and corrected imaginary coherence from complex
% cross-spectra
%
% This projects the cross-spectrum onto the imaginary axis to mitigate the
% effect of volume-conducted dependencies. Assumes volume-conducted sources
% arrive at sensors at the same time, resulting in a cross-spectrum with
% phase angle of 0 (perfectly in-phase) or pi (anti-phase) if the sensors
% are on opposite sides of a dipole source. With the imaginary coherence,
% in-phase and anti-phase associations are set to zero.
%
% The correction factor seeks to resolve artifically low connectivity near
% electrodes due to the imaginary component of the spectrum.
%
% Input:
% cs: either a matrix of size NxN for N sensors (or sources)
%     or a rank 3 or rank 4 tensor, where the first two indices
%     correspond to sensors. The third index should be frequency
%     and the fourth index could be time.
%
% Output:
% icoh: matrix (or matrices for tensors) of imaginary coherence with same
% dimsnsions as cs  
%
% cicoh: matrix (or matrices for tensors) of corrected imaginary coherence
% with same structure as cs
%
% Keeping in mind that the diagonal of the cross spectrum is the
% autospectrum, imaginary coherence is defined as:
%
%           Im(E[Sxy])
% C = ----------------------
%     sqrt(E[Sxx] * E[Syy])
%
% where Sxy is the cross spectrum, Sxx and Syy are the autospectra of
% channels X and Y, and E denotes the average across trials and Im denotes
% the imaginary part of the complex cross spectrum
%
% NOTE: ALL OF NOLTE'S WORK WAS BASED ON THE REAL AND IMAGINARY PARTS OF
% COHERENCY, AND NOT "COHERENCE" OR "MAGNITUDE SQUARED COHERENCE"
%
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Get dimensions and allocate memory for output:
%--------------------------------------------------------------------------
dims = size(cs);
ndim = numel(dims);
icoh = zeros(dims);
cicoh = zeros(dims);

switch ndim

    case 2
%--------------------------------------------------------------------------
% Channels x Channels (single frequency, no epochs):
%--------------------------------------------------------------------------  
        coy = cs./sqrt(diag(cs)*diag(cs)');
        %coyz = atan(coy);
        icoh = imag(coy);
        icoh(logical(eye(size(icoh)))) = 0;
        %icoh = atanh(icoh);
        rcoh = real(coy);
       
        for p = 1:dims(1)
           
            for q = 1:dims(1)
               
                if p ~= q
                    %cicoh(p,q) = icoh(p,q)/sqrt(1-rcoh(p,q).^2);
                    cicoh(p,q) = icoh(p,q)/sqrt(1-(rcoh(p,q)).^2);
                   
                end
           
            end
           
        end
       
    case 3
%--------------------------------------------------------------------------
% Channels x Channels x Frequency [or epochs] - depends on setup of cs:
%--------------------------------------------------------------------------
        [~, ~, nf] = size(cs);
        cicoh = zeros(62,62,nf);
       
        for f = 1:nf
            csloc = cs(:,:,f);
            coy = csloc./sqrt(diag(csloc)*diag(csloc)');
            %coyz = atanh(coy);
            ticoh = imag(coy);
            ticoh(logical(eye(size(ticoh)))) = 0;
            %ticoh = atan(ticoh);
            rcoh(:, :, f) = real(coy);
            icoh(:,:,f) = ticoh;
           
            for p = 1:dims(1)

                for q = 1:dims(1)

                    if p ~= q
                        cicoh(p,q,f) = icoh(p,q,f)/sqrt(1-rcoh(p,q,f).^2);

                    end

                end

            end
        end
       
    case 4
%--------------------------------------------------------------------------
% Channels x Channels x Frequency x Epochs:
%--------------------------------------------------------------------------
        [~, ~, nf, nep] = size(cs);
       
        for f = 1:nf
           
            for t = 1:nep
                csloc = cs(:,:,f,t);
                coy = csloc./sqrt(diag(csloc)*diag(csloc)');
                %coyz = atanh(coy);
                ticoh = imag(coy);
                ticoh(logical(eye(size(ticoh)))) = 0;
                %ticoh = atan(ticoh);
                icoh(:,:,f,t) = ticoh;
                rcoh(:, :, f, t) = real(coy);
               
                for p = 1:dims(1)
                   
                    for q = 1:dims(1)
                       
                        if p ~= q
                            cicoh(p,q,f,t) = icoh(p,q,f,t)/sqrt(1-rcoh(p,q,f,t).^2);
                       
                        end
           
                    end
           
                end
           
            end
           
        end
       
end

%--------------------------------------------------------------------------
% End function:
%--------------------------------------------------------------------------
end