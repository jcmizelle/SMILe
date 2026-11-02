function [mst_adj, node_dist, cost_value] = SMILe_MST_Euclidian(adj, thr)
%--------------------------------------------------------------------------
%
% This function computes the minimum spanning tree form of a graph through 
% the estimation of Euclidean distances (default). Other options include
% collelation and absloute value correlation.
%
% MST: use adj(n x p) to form (n-1) lines to connect (n) objects in the 
% shortest possible way in the (p) dimensional variable-space, under the 
% condition of 'no closed loops.'
% 
% Input:
% adj:          A nchan x nchan adjacency matrix from connectivity 
%               analysis. Can be weighted (non-binary) or binary. Must be 
%               undirected.
%
% thr:          Optional threshold. If not supplied, 95% CI will be used.
%               If the adjacency matrix has been previously thresholded, 
%               use 0 as input for this option.
%
% Output: 
% mst_adj:      Minimum spanning tree graph for the supplied adjacency
%               matrix. Will be of the same dimensions as adj.
%
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Check to see if a threshold (thr) was passed in. If not, calculate and 
% apply 2 standard deviations away from the mean as the threshold:
%-------------------------------------------------------------------------- 
if ~exist('thr', 'var') || isempty(thr)
    [CI, sd] = SMILe_CIcalc(adj);
    thr = sd(2,2);
    %thr = CI(2);
    adj(adj < thr) = 0;

%--------------------------------------------------------------------------
% Otherwise, use the user input:
%--------------------------------------------------------------------------
elseif exist('thr', 'var') && ~isempty(thr)
    adj(adj < thr) = 0;

end

%--------------------------------------------------------------------------
% Set mode:
%--------------------------------------------------------------------------
working_mode = 'euclidean';

%--------------------------------------------------------------------------
% Set exclusions, if any:
%--------------------------------------------------------------------------
exclude_adj = [];

%--------------------------------------------------------------------------
% Call the MST:
%--------------------------------------------------------------------------
[mst_adj,node_dist, cost_value] = mst_euclidean(adj,working_mode, exclude_adj);

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
function [mst_adj_out,node_dist_out, cost_value_out] = mst_euclidean(local_adj, working_mode, exclude_adj)

cost_value_out = 0;

if ~exist('working_mode')
    working_mode = 'euclidean';

end

if isempty(intersect({'euclidean','corr','abs_corr'}, working_mode))
    working_mode = 'euclidean';

end

if ~isempty(intersect({'corr','abs_corr'}, working_mode))
    local_adj = per_gene_normalization(local_adj); % this makes computing correlation easier
    local_adj = local_adj./norm(local_adj(1,:));

end

if ~exist('exclude_adj') || isempty(exclude_adj) 
    exclude_adj = sparse(size(local_adj,1),size(local_adj,1));

end

[nr, ~] = size(local_adj);

components = []; active_components =[];
mst_adj_out = sparse(size(local_adj,1),size(local_adj,1)); 
node_dist_out = sparse(size(local_adj,1),size(local_adj,1)); 
count = 0; 
fprintf('constructing a total of %d MST edges ... %6d', size(local_adj,1)-1,count);

for i = 1:nr
    if isequal(working_mode, 'euclidean')
        dist = comp_dist_euclidean(local_adj,i,1:nr); 
        
    elseif isequal(working_mode, 'corr')
        dist = comp_dist_corr(local_adj,i,1:nr);
        
    elseif isequal(working_mode, 'abs_corr')
        dist = comp_dist_abs_corr(local_adj,i,1:nr); 
        
    end
    
    dist(i) = max(dist)+1;
    dist = dist + exclude_adj(i,:).*(max(dist)+1);    
    dist = full(dist);
    [Dmin,Dwin] = min(dist);
    Xmst(i,:) = [i Dwin];
    
    if mst_adj_out(i,Dwin) ==0 && mst_adj_out(Dwin,i) == 0
        mst_adj_out(i,Dwin) = 1;
        mst_adj_out(Dwin,i) = 1;
        node_dist_out(i,Dwin) = Dmin;
        node_dist_out(Dwin,i) = Dmin;
        cost_value_out = cost_value_out + Dmin;
        count = count + 1; 
        fprintf('\b\b\b\b\b\b%6d', count);
    
    end
    
    if isempty(components)
        components = sparse(zeros(size(local_adj,1),1)); 
        components([i,Dwin],1) = 1; active_components = 1;
    
    else
        
        [existing_comp1] = find(components(i,:) == 1 & active_components == 1);
        [existing_comp2] = find(components(Dwin,:) == 1 & active_components == 1);
        
        if isempty(existing_comp1) && isempty(existing_comp2)
            components = [components,zeros(size(local_adj,1),1)]; 
            components([i,Dwin],end) = 1; 
            active_components = [active_components,1];
        
        elseif ~isempty(existing_comp1) && isempty(existing_comp2)
            components([i,Dwin],existing_comp1) = 1;
        
        elseif isempty(existing_comp1) && ~isempty(existing_comp2)
            components([i,Dwin],existing_comp2) = 1;
        
        elseif ~isempty(existing_comp1) && ~isempty(existing_comp2) && existing_comp1 ~= existing_comp2
            components = [components, components(:,existing_comp1)+components(:,existing_comp2)];
            active_components = [active_components,1];
            active_components([existing_comp1, existing_comp2]) = 0;
        
        end
        
    end
    
end
 
while sum(active_components)>1
    components_sizes = sum(components); components_sizes(active_components==0) = max(components_sizes+1);
    [dummy, existing_comp1] = min(components_sizes);
    ind1 = find(components(:,existing_comp1)==1); 
    ind1 = ind1(:)';
    ind2 = setdiff(1:size(components,1),ind1); 
    ind2 = ind2(:)';
    
    if isequal(working_mode, 'euclidean')
        dist = comp_dist_euclidean(local_adj,ind1,ind2); 
    
    elseif isequal(working_mode, 'corr')
        dist = comp_dist_corr(local_adj,ind1,ind2);
    
    elseif isequal(working_mode, 'abs_corr')
        dist = comp_dist_abs_corr(local_adj,ind1,ind2);
    
    end
    
    dist = dist + exclude_adj(ind1,ind2).*(max(max(dist))+1); dist = full(dist);
    [Dmin,ind] = min(reshape(dist,length(ind1)*length(ind2),1));
    j = ceil(ind/length(ind1));
    i = ind - (j-1)*length(ind1);
    Xmst = [Xmst; [ind1(i),ind2(j)]];
    mst_adj_out(ind1(i),ind2(j)) = 1; 
    mst_adj_out(ind2(j),ind1(i)) = 1; 
    node_dist_out(ind1(i),ind2(j)) = Dmin; 
    node_dist_out(ind2(j),ind1(i)) = Dmin; 
    cost_value_out = cost_value_out + Dmin;
    [existing_comp2] = find(components(ind2(j),:) == 1 & active_components == 1);
    components(:,existing_comp1) = components(:,existing_comp1) + components(:,existing_comp2);
    active_components(existing_comp2) = 0;
    count = count + 1; 
    fprintf('\b\b\b\b\b\b%6d', count);

end
fprintf('\n');

end

function dist = comp_dist_euclidean(X,ind1,ind2)
dist = zeros(length(ind1),length(ind2));
    
    for i=1:length(ind1)
        dist(i,:) = sqrt(sum((repmat(X(ind1(i),:),length(ind2),1) - X(ind2,:)).^2,2)); 
        
    end
    
end


function dist = comp_dist_corr(X,ind1,ind2)
    dist = zeros(length(ind1),length(ind2));
    corr = X(ind1,:)*X(ind2,:)';
    dist = 1-corr; 

end


function dist = comp_dist_abs_corr(X,ind1,ind2)
    dist = zeros(length(ind1),length(ind2));
    corr = X(ind1,:)*X(ind2,:)';
    dist = 1-abs(corr); 

end

%--------------------------------------------------------------------------
% Done with subroutines:
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% End for function:
%--------------------------------------------------------------------------
end
