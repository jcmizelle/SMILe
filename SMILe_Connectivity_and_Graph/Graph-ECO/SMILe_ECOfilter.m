function [binary, adjin] = SMILe_ECOfilter(adjin, directed)
%%
%--------------------------------------------------------------------------
%
% Efficiency Cost Optimization Filter
% 
% Threshold using efficiency cost optimization (ECO), to select a threshold 
% based on the optimization of the trade-off between the efficiency of a 
% network and its wiring cost.
%
% Reference: "A topological criterion to filter information in complex 
% brain networks, De Vico Fallani et al, Plos Comp Biol, 2017
%
% Input:
% - adj: weighted (directed or undirected) connectivity matrix
% - directed: is the matrix directed (1) on undirected (0)
% 
% Output:
% - binary: thresholded connectivity matrix
% - adj = keeps the weights that overcome the threshold   
%
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Make a copy of the original matrix to retain Wrights:
%--------------------------------------------------------------------------
adj = adjin;

N = size(adj,1);

%--------------------------------------------------------------------------
% Determine if directed, and reshape if so:
%--------------------------------------------------------------------------
if directed
    numcon = 3*N;
    ind = find(adj);
    
else
    adj = triu(adj);
    numcon = 1.5*N;
    ind = find(triu(adj));
    
end

%--------------------------------------------------------------------------
% Check quality of the matrix:
%--------------------------------------------------------------------------
if numcon>length(ind)
    error('Input matrix is too sparse')
    
end

%--------------------------------------------------------------------------
% Start ECO:
%--------------------------------------------------------------------------
sorind = sortrows([ind adj(ind)],-2);

adj(sorind(numcon+1:end,1)) = 0;

if directed
    binary = double(logical(adj));
    
else
    binary = double(logical(adj+adj'));
    
end

%--------------------------------------------------------------------------
% Apply threshold to weighted matrix:
%--------------------------------------------------------------------------
adjin(binary ~= 1) = 0;

%--------------------------------------------------------------------------
% End for function:
%--------------------------------------------------------------------------
end