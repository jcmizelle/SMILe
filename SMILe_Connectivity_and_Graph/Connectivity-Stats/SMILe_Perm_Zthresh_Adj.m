function adjthr = SMILe_Perm_Zthresh_Adj(adj, cfg)
%--------------------------------------------------------------------------
%
% This script generates a thresholded adjacency matrix based on permutation
% statistics and z-statistic calculations. 
%
% Input:
% adj:      An adjacency matrix to permute and threshold
%
% cfg:      A configuration structure with options
% nperms:   Number of permutations. Default is 1000.
% pval:     P-value for thresholding the adjacency matrix. Default is 0.05.
% doplot:   Option to generate our line plot with the thresholded adjacency
%               matrix. Default is 0 (no plot).
%
% Output:
% adjthr:   nchan x nchan adjacaney matrix thresholded at the requested
%           level.
%01
% Example:
% cfg = [];
% cfg.nperms = 1000;
% cfg.pval = 0.05;
% cfg.doplot = 0;
% cfg.layout = layout;
% adjthr = SMILe_Perm_Zthresh_Adj(adj, cfg)
%
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Testing mechanics:
%--------------------------------------------------------------------------
% adjmin = -1;
% adjmax = 1;
% adj = (adjmax-adjmin).*rand(62, 62) + adjmin;

%--------------------------------------------------------------------------
% Parse input:
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Set defaults number of permutations if not present in cfg:
%--------------------------------------------------------------------------
if isfield(cfg, 'nperms')
    nperms = cfg.nperms;
    
else
    nperms = 1000;

end

%--------------------------------------------------------------------------
% Set p-value level if not in cfg (0 = no, 1 = yes):
%--------------------------------------------------------------------------
if isfield(cfg, 'pval')
    pval = cfg.pval;
    
else
    pval = 0.05;

end

%--------------------------------------------------------------------------
% Set plot option if not in cfg (0 = no, 1 = yes):
%--------------------------------------------------------------------------
if isfield(cfg, 'doplot')
    doplot = cfg.doplot;
    
else
    doplot = 0;

end

%--------------------------------------------------------------------------
% Get dimensions of the adjacency matrix:
%--------------------------------------------------------------------------
dims = size(adj);

%--------------------------------------------------------------------------
% Take care of adj properties (get absolute value, make symmetrical, set 
% diagonals to zero):
%--------------------------------------------------------------------------
adj = abs(adj);
adj = triu(adj);
adj = adj + adj';
adj(logical(eye(size(adj)))) = 0;

%--------------------------------------------------------------------------
% Allocate memory to store permutations of adj:
%--------------------------------------------------------------------------
permdat = zeros(dims(1), dims(2), nperms);

%--------------------------------------------------------------------------
% Start permutation process:
%--------------------------------------------------------------------------
for np = 1:nperms
    permidxr = randperm(dims(1));
    permidxc = randperm(dims(2));
    permdat(:,:,np) = adj(permidxr, permidxc);
    
end

%--------------------------------------------------------------------------
% Compute z-map:
%--------------------------------------------------------------------------
adjz = (adj-mean(permdat, 3))./std(permdat, [], 3);

%--------------------------------------------------------------------------
% Apply z-map threshold to the original adjacency matrix:
%--------------------------------------------------------------------------
adjthr = adj;
adjthr(abs(adjz) < norminv(1-pval)) = 0;

%--------------------------------------------------------------------------
% Plot if requested:
%--------------------------------------------------------------------------
if doplot == 1 && isfield(cfg, 'layout') && ~isempty(cfg.layout)

    try 
        layout = cfg.layout;
        SMILe_Line_Plot(adjthr, layout)

    catch
        %disp('No Layout')
        SMILe_Line_Plot(adjthr, [])
        
    end
    

end

%--------------------------------------------------------------------------
% End function:
%--------------------------------------------------------------------------
end