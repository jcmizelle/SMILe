function [binary, thr_adj] = SMILe_MST_Prim(adj)
%--------------------------------------------------------------------------
%
% This function computes the minimum spanning tree form of a graph through 
% the estimation of Prim's algorithm.
%
% Prim's basic idea: start at any node, find closest neighbor and mark 
% edges for all remaining nodes, find closest to previous cluster, mark 
% edge continue until no nodes remain.
% 
% Input:
% adj:          A nchan x nchan adjacency matrix from connectivity 
%               analysis. Must be weighted (non-binary). Must be 
%               undirected.
%
% Output: 
% mst_adj:      Minimum spanning tree graph for the supplied adjacency
%               matrix. Will be of the same dimensions as adj.
%
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Ensure symmetry and connectedness:
%--------------------------------------------------------------------------
adj = triu(adj);
adj = adj+adj';

%--------------------------------------------------------------------------
% Call the MST:
%--------------------------------------------------------------------------
binary = mst_Prim(adj);

%--------------------------------------------------------------------------
% Apply threshold to weighted matrix:
%--------------------------------------------------------------------------
thr_adj = adj;
thr_adj(binary ~= 1) = 0;


%--------------------------------------------------------------------------
% Create graph object:
%--------------------------------------------------------------------------
% nodelocs  = elec.chanpos(:,1:2);
% nodenames = elec.label;
% 
% G = [];
% G = graph(mst_adj, 'upper', 'omitselfloops');
% if size(mst_adj,1) < size(adj,1)
%     G = addnode(G,size(adj,1) - size(mst_adj,1));
%     
% end
% 
% G.Nodes.Name = nodenames';
% G.Nodes.Loc = nodelocs; 

%--------------------------------------------------------------------------
% Begin subroutines:
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Start mst_Prim:
%--------------------------------------------------------------------------
function mst_adj_out = mst_Prim(adj)

% check if graph is connected:
if ~(local_isconnected(adj))
    fprintf('the graph is not connected, no spanning tree exists\n'); 
    
    return 

end

n = length(adj); % number of nodes
mst_adj_out = zeros(n);   % initialize tree

adj(find(adj == 0)) = inf; % set all zeros in the matrix to inf

conn_nodes = 1;        % nodes part of the min-span-tree
rem_nodes = [2:n];     % remaining nodes

while length(rem_nodes) > 0
  [minlink] = min(min(adj(conn_nodes,rem_nodes)));
  ind = find(adj(conn_nodes,rem_nodes) == minlink);

  [ind_i,ind_j] = ind2sub([length(conn_nodes),length(rem_nodes)],ind(1));

  i = conn_nodes(ind_i); 
  j = rem_nodes(ind_j); % gets back to adj indices
  mst_adj_out(i,j) = 1; 
  mst_adj_out(j,i) = 1;
  conn_nodes = [conn_nodes j];
  rem_nodes = setdiff(rem_nodes, j);
  
end

%--------------------------------------------------------------------------
% End MST:
%--------------------------------------------------------------------------
end

%--------------------------------------------------------------------------
% Start isconnected:
%--------------------------------------------------------------------------
function S = local_isconnected(adj)

if ~(isempty(find(sum(adj) == 0))) 
    S = false; 
    return 

end

n = length(adj);
x = [1; zeros(n-1,1)]; % [1,0,...0] nx1 vector 

while 1
     y = x;
     x = adj*x + x;
     x = x > 0;
     
     if x == y
         break 
     
     end

end

S = true;

if sum(x) < n
    S = false; 

end

%--------------------------------------------------------------------------
% End isconnected:
%--------------------------------------------------------------------------
end

%--------------------------------------------------------------------------
% End subroutines:
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% End function:
%--------------------------------------------------------------------------
end