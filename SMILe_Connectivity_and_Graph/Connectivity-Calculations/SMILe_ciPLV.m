function [ciPLV] = SMILe_ciPLV(cs)
%--------------------------------------------------------------------------
% The phase locking value attempts to mitigate power differences
% between realizations (tapers or trials) by treating all values of
% the cross-spectrum as the same power. This has the effect of
% downweighting high power realizations and upweighting low power
% realizations.
%
% Input: 
% cs: The (complex valued) cross-spectrum for in the format of an
% nchan x nchan x epochs matrix 
%
% PLV = |E[Sxy/|Sxy|]|
%
% iPLV = |E[Im(Sxy/|Sxy|)]|
%
%                  |E[Im(Sxy/|Sxy|)]|
% ciPLV = ------------------------------------
%          sqrt(1 - |E[real(Sxy/|Sxy|)]| ^ 2)
%
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Imaginary and Real PLV:
%--------------------------------------------------------------------------
for nep = 1:size(cs,3)
    epiPLV(:,:,nep) = imag(cs(:,:,nep)./abs(cs(:,:,nep)));
    eprPLV(:,:,nep) = real(cs(:,:,nep)./abs(cs(:,:,nep)));
    
end

%--------------------------------------------------------------------------
% Mean Imaginary PLV:
%--------------------------------------------------------------------------
siPLV = abs(sum(epiPLV,3));
miPLV = siPLV / nep;

%--------------------------------------------------------------------------
% Mean Real PLV:
%--------------------------------------------------------------------------
srPLV = abs(sum(eprPLV,3));
mrPLV = srPLV / nep;

%--------------------------------------------------------------------------
% Corrected Imaginary PLV:
%--------------------------------------------------------------------------
ciPLV =  miPLV ./ sqrt(1 - (mrPLV.^2));

%--------------------------------------------------------------------------
% Force symmetry:
%--------------------------------------------------------------------------
ciPLV = triu(ciPLV);
ciPLV = ciPLV + ciPLV';

%--------------------------------------------------------------------------
% Set diagonals to zero:
%--------------------------------------------------------------------------
ciPLV(logical(eye(size(ciPLV)))) = 0;

%--------------------------------------------------------------------------
% End function:
%--------------------------------------------------------------------------
end