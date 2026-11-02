function [pPLV] = SMILe_partialPLV(cs)
%--------------------------------------------------------------------------
%
% Calculates partial phase locking value from cross-spectra
%
% This essentially represents the inversion of a matrix of pairwise PLVs to 
% compute a partial phase locking value. 
% 
% This is analogous to inversion of cross-correlation and cross-spectral 
% matrices to compute partial correlations and partial coherence.
%
% B. Schelter et al. (2006), Partial phase synchronization for 
%   multivariate synchronizing system, Phys. Rev. Lett., 96, 208103. 
%
% Sergul Aydore, Dimitrios Pantazis and Richard Leahy, Measuring 
%   Partial Phase Locking Value to Detect Synchronization in Multivariate 
%   Gaussian Systems
% 
% Input:
% cs:       A complex valued cross spectrum matrix of size nchan x nchan.
%
%           cs is assumed to be of a single frequency or an averaged band
%
% Output: 
% pPLV:     A matrix of partial phase locking values with same structure 
%           as cs 
%
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Check dimensions and allocate memory:
%--------------------------------------------------------------------------
dims = size(cs);
pPLV = zeros(dims);

%--------------------------------------------------------------------------
% Start calculations:
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Calculate the pseudoinverse:
%--------------------------------------------------------------------------
invConn = pinv(cs);

%--------------------------------------------------------------------------
% Grab the diagonals and calculate the delimiter:
%--------------------------------------------------------------------------
autoInv = diag(invConn);
delim = autoInv * autoInv';

%--------------------------------------------------------------------------
% Compute the partial phase locking value:
%--------------------------------------------------------------------------
pPLV = real(abs(invConn) ./ sqrt(delim));

%--------------------------------------------------------------------------
% Set diagonals to zero:
%--------------------------------------------------------------------------
pPLV(logical(eye(size(pPLV)))) = 0;

end