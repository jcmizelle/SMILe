function [adjthr1, adjthr2] = SMILe_Perm_Zthresh_DiffAdj(adj1, adj2, cfg)
%--------------------------------------------------------------------------
%
% This script generates a thresholded adjacency matrix based on permutation
% statistics and z-statistic calculations. 
%
% Input:
% adj1:     An array of adjacency matrices (nchan x nchan x nsub) for
%           condition or group 1
% adj2:     An array of adjacency matrices (nchan x nchan x nsub) for
%           condition or group 2
%
% cfg:      A configuration structure with options
%   nperms:   Number of permutations. Default is 1000.
%   pval:     P-value for thresholding the adjacency matrix. Default is 
%             0.05.
%
% Output:
% adjthr1:  nchan x nchan adjacaney matrix thresholded at the requested
%           level. This should relate to adj1 > adj2.
% adjthr2:  nchan x nchan adjacaney matrix thresholded at the requested
%           level. This should relate to adj2 > adj1.
%
% Example:
% cfg = [];
% cfg.nperms = 1000;
% cfg.pval = 0.05;
% adjthr = SMILe_Perm_Zthresh_Adj(adj, cfg)
%
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Testing mechanics:
%--------------------------------------------------------------------------
% adjmin = eps;
% adjmax = .5;
% adj1 = (adjmax-adjmin).*rand(62, 62, 15) + adjmin;
% adj2 = (adjmax-adjmin).*rand(62, 62, 15) + adjmin;

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
dimsadj1 = size(adj1);
dimsadj2 = size(adj2);

%--------------------------------------------------------------------------
% Take care of adj properties (get absolute value, make symmetrical, set 
% diagonals to zero):
%--------------------------------------------------------------------------
for d = 1:dimsadj1(end)
    tmpadj = squeeze(adj1(:,:,d));
    tmpadj = abs(tmpadj);
    tmpadj = tmpadj + tmpadj';
    tmpadj(logical(eye(size(tmpadj)))) = 0;
    adj1(:,:,d) = tmpadj;
    
end

clear tmpadj

for d = 1:dimsadj2(end)
    tmpadj = squeeze(adj2(:,:,d));
    tmpadj = abs(tmpadj);
    tmpadj = tmpadj + tmpadj';
    tmpadj(logical(eye(size(tmpadj)))) = 0;
    adj2(:,:,d) = tmpadj;
    
end

clear tmpadj

%--------------------------------------------------------------------------
% Concatenate groups/conditions/etc:
%--------------------------------------------------------------------------
adj = cat(3, adj1, adj2);

%--------------------------------------------------------------------------
% Get dimensions of concatenated data:
%--------------------------------------------------------------------------
dims = size(adj);

%--------------------------------------------------------------------------
% Start permutation process. If even groups, calculate actual differences
% of individual permutations:
%--------------------------------------------------------------------------
if dimsadj1(end) == dimsadj2(end)

%--------------------------------------------------------------------------
% Get true difference (of individuals):
%--------------------------------------------------------------------------
truediff = adj1 - adj2;

%--------------------------------------------------------------------------
% Allocate memory to store permutations of adj:
%--------------------------------------------------------------------------
permdiff = zeros(dims(1), dims(2), nperms, dims(3)*.5);

%--------------------------------------------------------------------------
% Start permutations:
%--------------------------------------------------------------------------
    for np = 1:nperms
        permidxsub = randperm(dims(3));
        permdiff(:,:, np, :) = adj(:, :, permidxsub(1:dims(3)*.5)) - adj(:, :, permidxsub(dims(3)*.5+1:end));
        
    end

%--------------------------------------------------------------------------
% Take average of permuted differences:
%--------------------------------------------------------------------------
    permdiff = mean(permdiff, 4);

%--------------------------------------------------------------------------
% To totally scramble the data:
%--------------------------------------------------------------------------
%     for np = 1:nperms
%         permidxr = randperm(dims(1));
%         permidxc = randperm(dims(2));
%         permsub = randperm(dims(3));
%         permdiff(:,:,np) = mean(adj(permidxr, permidxc, permsub(1:dims(3)*.5)), 3) - mean(adj(permidxr, permidxc, permsub(dims(3)*.5+1:end)), 3);
% 
%     end

%--------------------------------------------------------------------------
% If uneven groups, calculate permuted group mean differences:
%--------------------------------------------------------------------------
elseif dimsadj1(end) ~= dimsadj2(end)

%--------------------------------------------------------------------------
% Get true difference (of means):
%--------------------------------------------------------------------------
    truediff = mean(adj1, 3) - mean(adj2, 3);

%--------------------------------------------------------------------------
% Start permutations:
%--------------------------------------------------------------------------
    for np = 1:nperms
        permidxr = randperm(dims(1));
        permidxc = randperm(dims(2));
        permidxsub = randperm(dims(3));
        permdiff(:,:,np) = mean(adj(permidxr, permidxc, permidxsub(1:dimsadj1(end))), 3) - mean(adj(permidxr, permidxc, permidxsub(dimsadj1(end)+1:end)), 3);

    end
    
end

%--------------------------------------------------------------------------
% Compute z-map:
%--------------------------------------------------------------------------
adjz = (truediff-mean(permdiff, 3))./std(permdiff, [], 3);

%--------------------------------------------------------------------------
% Apply z-map threshold to the original adjacency matrix:
%--------------------------------------------------------------------------
adjthr = truediff;
adjthr(abs(adjz) < norminv(1-pval)) = 0;

%--------------------------------------------------------------------------
% Take direction of difference into account:
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% For Cond 1 > Cond 2:
%--------------------------------------------------------------------------
adjthr1(adjthr < 0) = 0;

%--------------------------------------------------------------------------
% For Cond 2 > Cond 1:
%--------------------------------------------------------------------------
adjthr2(adjthr > 0) = 0;

%--------------------------------------------------------------------------
% End function:
%--------------------------------------------------------------------------
end