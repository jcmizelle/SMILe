function [G] = SMILe_adj2graph(adj, elec)
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
% COHERENCY, AND NOT 
%
%--------------------------------------------------------------------------

nodelocs  = elec.chanpos;
nodenames = elec.label;

%--------------------------------------------------------------------------
% Check whether the adjacency matrix is binary:
%--------------------------------------------------------------------------
[~, idx] = find(adj > 0);

%--------------------------------------------------------------------------
% adj cannot be binary if its sum isn't equal to the number of nonzero 
% elements:
%--------------------------------------------------------------------------
if sum(sum(adj)) ~= numel(idx)

%--------------------------------------------------------------------------
% Check to see if a threshold (thr) was passed in. If not, calculate and 
% apply 2 standard deviations away from the mean as the threshold:
%-------------------------------------------------------------------------- 
%     if ~exist('thr', 'var') || isempty(thr)
%         [~, sd] = SMILe_CIcalc(adj);
%         thr = sd(2,2);
%         %thr = CI(2);
%         adj(adj >= thr) = 1;
%         adj(adj < thr) = 0;
% 
% %--------------------------------------------------------------------------
% % Otherwise, use the user input:
% %--------------------------------------------------------------------------     
%     else
%         adj(adj >= thr) = 1;
%         adj(adj < thr) = 0;
% 
%     end

end

%--------------------------------------------------------------------------
% Build the graph object:
%--------------------------------------------------------------------------
G = [];
G = graph(adj, nodenames, 'upper', 'omitselfloops');
G.Nodes.Name = nodenames';
G.Nodes.Loc = nodelocs; 

%--------------------------------------------------------------------------
% End for function:
%--------------------------------------------------------------------------
end