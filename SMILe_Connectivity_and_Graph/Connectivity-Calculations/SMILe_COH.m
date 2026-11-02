function [mscoh, coh, coy] = SMILe_COH(cs)
%--------------------------------------------------------------------------
%
% Calculates coherency and magnitude squared coherency from complex 
% cross-spectra. Coherence represents the complex-valued linear association 
% between time series in the frequency domain.
%
% Coherency is the normalized cross-spectrum at a specific frequency, i.e. 
% for two signals A and B it is defined as Cab/sqrt(Caa * Cbb). Coherency 
% is a complex-valued number with a magnitude (or modulus) smaller than 
% or equal to 1 and an angle corresponding to the phase lag between the 
% signals.
%
% The modulus of the coherency is known as the coherence. The coherence 
% is a real-valued number between 0 and 1. Note that it is also not 
% uncommon in scientific publications to report the magnitude-squared 
% coherence.
% 
% Input:
% cs:       either a matrix of size NxN for N sensors (or sources) 
%           or a rank 3 or rank 4 tensor, where the first two indices
%           correspond to sensors. The third index should be frequency 
%           and the fourth index could be time/epochs.
%
% Output: 
% coy:      matrix (or matrices for tensors) of complex coherency with same 
%           structure as cs  
% 
% coh:      matrix (or matrices for tensors) of the real part of complex 
%           coherency with same structure as cs  
%
% mscoh:    matrix (or matrices for tensors) of complex magnitude squared 
%           coherence with same structure as cs
%
%
% Keeping in mind that the diagonal of the cross spectrum is the
% autospectrum, coherence is defined as:
%
%         {      | E[Sxy] | ^2    }
% MSCOH = { --------------------- }
%         { sqrt(E[Sxx] * E[Syy]) }
%
%              E[Sxy]     
% COH = ---------------------
%        sqrt(E[Sxx] * E[Syy])
%
% where Sxy is the cross spectrum, Sxx and Syy are the autospectra of
% channels X and Y, and E denotes the average across trials
%
%--------------------------------------------------------------------------


%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Get started:
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

%--------------------------------------------------------------------------
% Get dimensions and allocate memory for output:
%--------------------------------------------------------------------------
dims = size(cs);
ndim = numel(dims);
mscoh = zeros(dims);
coh = zeros(dims);
coy = zeros(dims);

%--------------------------------------------------------------------------
% Build conditional on dimensions:
%--------------------------------------------------------------------------
switch ndim
    
    case 2
%--------------------------------------------------------------------------
% Channels x Channels (single frequency, no epochs):
%--------------------------------------------------------------------------        
        mscoh = (abs(cs).^2)./(diag(cs)*diag(cs)');
        mscoh(logical(eye(size(mscoh)))) = 0;
        
        coy = cs./sqrt(diag(cs)*diag(cs)');
        coy(logical(eye(size(coy)))) = 0;    
        
    case 3
%--------------------------------------------------------------------------
% Channels x Channels x Frequency [or epochs] - depends on setup of cs:
%--------------------------------------------------------------------------        
        [nc, nc, nf] = size(cs);
        coh = cs;
        
        for f = 1:nf
            csloc = squeeze(cs(:,:,f));
            
            tmscoh = (abs(csloc).^2)./(diag(csloc)*diag(csloc)');
            tmscoh(logical(eye(size(tmscoh)))) = 0;
            mscoh(:,:,f) = mstcoh;
            
            tcoy = csloc./sqrt(diag(csloc)*diag(csloc)');
            %tcoh = tcoh;
            tcoy(logical(eye(size(tcoy)))) = 0;
            coy(:,:,f) = tcoy;
        
        end
        
    case 4
%--------------------------------------------------------------------------
% Channels x Channels x Frequency x Epochs:
%--------------------------------------------------------------------------
        [nc, nc, nf, nt] = size(cs);
        coh = cs;
        
        for f = 1:nf
            
            for t = 1:nt
                csloc = squeeze(squeeze(cs(:,:,f,t)));
                
                tmscoh = (abs(csloc).^2)./(diag(csloc)*diag(csloc)');
                tmscoh(logical(eye(size(tmscoh)))) = 0;
                mscoh(:,:,f,t) = tmscoh;
                
                tcoy = abs(csloc)./sqrt(diag(csloc)*diag(csloc)');
                %tcoh = tcoh;
                tcoy(logical(eye(size(tcoy)))) = 0;
                coy(:,:,f,t) = tcoy;
            
            end
            
        end
        
end

%--------------------------------------------------------------------------
% Take the real part of the complex cross spectrum for real coherence and 
% magnitude squared coherence:
%-------------------------------------------------------------------------- 
mscoh = real(mscoh);
coh = real(coy);

%--------------------------------------------------------------------------
% End function:
%-------------------------------------------------------------------------- 
end