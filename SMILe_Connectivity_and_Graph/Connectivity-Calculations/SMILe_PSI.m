function [psi, psi_z] = SMILe_PSI(cs, freqs)
%--------------------------------------------------------------------------
%
% This function alculates the phase slope index (PSI), a measurement of
% effective connectivity.
%
% Usage:
%   [psi, psi_z] = SMILe_PSI(cs, freqbins)
%
% Input:
% cs:           Complex valued cross spectrum calculated either through
%               SMILe_CS_Fourier or SMILe_CS_Wavelet.
% freqs:        M x 2 matrix containing frequency boundaries in rows M.
%               Example: [4 8; 8 12; 13 30]
% permtest:     Permutation test (boolean). If true, psi output values are
%               in standardized Z values relative to a null hypothesis
%               distribution
%
% Output:
% psi:          phase slope index values. For M channels PSI is an M x M
%               matrix if one frequency bin, or M x M x K if freqs has K
%               rows (with K>1).
%              
%               psi(i,j) is the directed connectivity from channel i to
%               channel j, (e.g., channel i is the sender if psi(i,j) is
%               positive; channel i is the receiver if negative)
%
% Computed through methods formulated in:
% Nolte G, Ziehe A, Nikulin VV, Schlogl A, Kramer N, Brismar T, Muller KR.
% Robustly estimating the flow direction of information in complex physical
% systems. Physical Review Letters.
% See also: http://doc.ml.tu-berlin.de/causality/
%
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Begin PSI calculations
%--------------------------------------------------------------------------
freqbins = freqs + 1;
nchan = size(cs,1);
csall = cs;

%for freqbini=1:size(freqbins,1)
   
%--------------------------------------------------------------------------
% Find FFT indices of requested frequency bands:
%--------------------------------------------------------------------------
[nchan nchan nf] = size(cs);

[nm nf] = size(freqbins);
psi = zeros(nchan, nchan, nm);
psisum = zeros(nchan, nm);

for ii = 1:nm;
    psi(:,:,ii) = cs2ps(cs(:,:,freqbins(ii,:)));
    psisum(:,ii) = sum(psi(:,:,ii),2);

end
    
%--------------------------------------------------------------------------
% Set diagonals to zero:
%--------------------------------------------------------------------------
psi(logical(eye(size(psi)))) = 0;
psisum(logical(eye(size(psisum)))) = 0;
   
   
%--------------------------------------------------------------------------
% Mandatory permutation testing with 10,000 permutations:
%--------------------------------------------------------------------------
n_permutes = 1000;

%--------------------------------------------------------------------------
% Create null distribution:
%--------------------------------------------------------------------------
nulldist = zeros(nchan,nchan,n_permutes);
nfidx = randi(size(cs,3),4,1);

for freqbini = 1:nm
    for permi = 1:n_permutes
        nulldist(:,:,permi) = cs2ps(cs(:,:,nfidx));
        %nulldist(:,:,permi) = sum(imag(conj(psi(:,:,randperm(nfidx))).*psi(:,:,randperm(nfidx))),3);

    end

psi_z(:,:,freqbini) = (psi(:,:,freqbini) - mean(nulldist,3)) ./ std(nulldist,[],3);

%--------------------------------------------------------------------------
% Set diagonals to zero:
%--------------------------------------------------------------------------
tmp = squeeze(psi_z(:,:,freqbini));
tmp(logical(eye(size(tmp)))) = 0;
psi_z(:,:,freqbini) = tmp;

end
%
%--------------------------------------------------------------------------
% Subroutine to create the PSI :
%--------------------------------------------------------------------------
function ps = cs2ps(cs);
        df = 1;
        [nchan nchan nf]=size(cs);
        pp = cs;
        for f = 1:nf
            pp(:,:,f) = cs(:,:,f)./sqrt(diag(cs(:,:,f))*diag(cs(:,:,f))');
        
        end
        
        ps = sum(imag(conj(pp(:,:,1:end-df)).*pp(:,:,1+df:end)),3);

end
%--------------------------------------------------------------------------
% End for function:
%--------------------------------------------------------------------------
end