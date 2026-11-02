function [PLV] = SMILe_PLV(cs)
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
%                  |E[Im(Sxy/|Sxy|)]|
% ciPLV = ------------------------------------
%          sqrt(1 - |E[real(Sxy/|Sxy|)]| ^ 2)
%
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Standard PLV:
%--------------------------------------------------------------------------
for nep = 1:size(cs,3)
    epPLV(:,:,nep) = cs(:,:,nep)./abs(cs(:,:,nep));
    
end

%--------------------------------------------------------------------------
% Get the mean of the PLV:
%--------------------------------------------------------------------------
sPLV = abs(sum(epPLV,3));
PLV = sPLV./nep;

%--------------------------------------------------------------------------
% Force symmetry:
%--------------------------------------------------------------------------
PLV = triu(PLV);
PLV = PLV + PLV';

%--------------------------------------------------------------------------
% Set diagonals to zero:
%--------------------------------------------------------------------------
PLV(logical(eye(size(PLV)))) = 0;

%--------------------------------------------------------------------------
% End function:
%--------------------------------------------------------------------------
end