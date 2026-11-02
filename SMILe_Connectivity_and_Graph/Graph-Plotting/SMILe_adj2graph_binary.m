function [G] = SMILe_adj2graph_binary(adj, elec, thr)
%--------------------------------------------------------------------------
%
% This function creates a graph object (binary) from an adjacency matrix
% that was calculated through any of the SMILe_XXX connectivity
% calculations.
% 
% Input:
% adj:  A nchan x nchan adjacency matrix from connectivity analysis. Can be 
%       weighted (non-binary) or binary. Must be undirected.
% 
% elec: The 64-channel electrode position template
%
% thr:  Optional threshold. If not supplied, +2 SD will be used
%
% Output: 
% G:    Graph object for display
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
    if ~exist('thr', 'var') || isempty(thr)
        [CI, sd] = SMILe_CIcalc(adj);
        thr = sd(2,2);
        %thr = CI(2);
        adj(adj >= thr) = 1;
        adj(adj < thr) = 0;

%--------------------------------------------------------------------------
% Otherwise, use the user input:
%--------------------------------------------------------------------------     
    else
        adj(adj >= thr) = 1;
        adj(adj < thr) = 0;

    end

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