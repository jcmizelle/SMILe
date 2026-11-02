function SMILe_Circuit_Plot(G, rminf)
%--------------------------------------------------------------------------
%
% This function creates a graph object (binary) from an adjacency matrix
% that was calculated through any of the SMILe_XXX connectivity
% calculations.
% 
% Input:
% G:        Graph object from SMILe_adj2graph_binary.
% 
% rminf:    Option to NOT display nodes with no connecting edges.
%
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Check if non-connected nodes should be included:
%--------------------------------------------------------------------------
if exist('rminf', 'var') && ~isempty(rminf) && rminf > 0
    G = rmnode(G, find(degree(G) == 0)');
    
end

%--------------------------------------------------------------------------
% Start figure, with the following layout options:
% 
% 'circle', 'force', 'layered', 'subspace', 'force3', or 'subspace3'
%--------------------------------------------------------------------------
figure;
fig = plot(G,'Layout','layered');

%--------------------------------------------------------------------------
% Force user options to figure:
%--------------------------------------------------------------------------
fig.Marker = 'o';
fig.NodeColor = "#C8A2C8";
fig.MarkerSize = 5;
fig.LineWidth = 1.5;
fig.EdgeColor = "#DCD0FF";
nl = fig.NodeLabel;
fig.NodeLabel = '';
nudgefactor = 0.5;
xd = get(fig, 'XData') + nudgefactor;
yd = get(fig, 'YData');
text(xd, yd, nl, 'FontSize',10, 'FontWeight','bold', 'HorizontalAlignment','left', 'VerticalAlignment','middle')

%--------------------------------------------------------------------------
% End for function:
%--------------------------------------------------------------------------
end