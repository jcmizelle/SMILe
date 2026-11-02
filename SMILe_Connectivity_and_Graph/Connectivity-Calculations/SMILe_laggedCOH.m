function [lagCOH] = SMILe_laggedCOH(cs)
%--------------------------------------------------------------------------
%
% Calculates lagged coherence from the complex cross-spectrum. 
%
% Pascual-Marqui et al. Assessing interactions in the brain with 
%   exact low-resolution electromagnetic tomography Phil. Trans. R. Soc. A 
%   2011 369, 3768-3784
% 
%
% Input: 
% cs:   nchan x nchan matrix of the complex valued cross-spectrum
%
%           cs is assumed to be of a single frequency or an averaged band
%
% Output
% lagCOH:  nchan x nchan matrix of lagged coherence
%
% These calculations are absed off of Guido Nolte's formulation of lagged
% coherence.
%
%                Im(E[Sxy])
% C = --------------------------------
%     sqrt(1 - real(E[Sxx] * E[Syy])^2 + eps)
%
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Check dimensions:
%--------------------------------------------------------------------------
dims = size(cs);
lagCOH = zeros(dims);

%--------------------------------------------------------------------------
% Compute coherency:
%--------------------------------------------------------------------------
coh = cs./sqrt(diag(cs)*diag(cs)');

%--------------------------------------------------------------------------
% Calculte lagged coherence:
%--------------------------------------------------------------------------
lagCOH = imag(coh)./sqrt((1-real(coh).^2+eps));

%--------------------------------------------------------------------------
% Set diagonals to zero:
%--------------------------------------------------------------------------
lagCOH(logical(eye(size(lagCOH)))) = 0;

%--------------------------------------------------------------------------
% End function:
%-------------------------------------------------------------------------- 
end