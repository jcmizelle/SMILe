function [iPLV] = SMILe_iPLV(cs)
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
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Get the imaginary part of the PLV:
%--------------------------------------------------------------------------
for nep = 1:size(cs,3)
    epiPLV(:,:,nep) = imag(cs(:,:,nep)./abs(cs(:,:,nep)));
    
end

%--------------------------------------------------------------------------
% Get the mean of the imaginary PLV:
%--------------------------------------------------------------------------
siPLV = abs(sum(epiPLV,3));
iPLV = siPLV / nep;

%--------------------------------------------------------------------------
% Force symmetry:
%--------------------------------------------------------------------------
iPLV = triu(iPLV);
iPLV = iPLV + iPLV';

%--------------------------------------------------------------------------
% Set diagonals to zero:
%--------------------------------------------------------------------------
iPLV(logical(eye(size(iPLV)))) = 0;

%--------------------------------------------------------------------------
% End function:
%--------------------------------------------------------------------------
end