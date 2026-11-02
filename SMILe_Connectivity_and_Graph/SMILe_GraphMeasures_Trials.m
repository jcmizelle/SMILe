%--------------------------------------------------------------------------        
% This is SMILe's code to calculate graph metrics. It will accept an
% adjacency matrix of varying sizes (see below), but the matrix should have
% come through SMILE's cross spectral calculations (either Fourier or
% Morelet Wavelet). Many of the function calls are based on the Brain
% Connectivity Toolbox, and others are custom.
%
% There are many metrics calculated here - see a description of each in the
% script SMILe_GraphMeasuresDescribed.
% 
% Also note that we are calculating both nodal and global measures, and so
% your eventual plotting and statistical processes will need to be adjusted
% accordingly.
%
% I have not run into this yet, but there is always the chance that you
% will need to "squeeze" the output of the function to take care of
% singleton dimensions.
%--------------------------------------------------------------------------
function GraphMeasuresOut = SMILe_GraphMeasures_Trials(mst_adj_in)

%--------------------------------------------------------------------------        
% Parse input for multiple dimensions. If multidimensional and coming from 
% SMILe's TF-Wavelet code that calculates things, then the data should be 
% in the format of: nchan x nchan x time x subjects

% 11/12/2021: Editing the file to allow nchan x nchan x time x subjects x trials
%
% TODO: build in an option for multiple bands
% TODO: build a cell structure for hub calculations which are always of 
% different dimensions
%--------------------------------------------------------------------------
mst_adj_in = full(mst_adj_in);
dims = size(mst_adj_in);

%--------------------------------------------------------------------------        
% Full (non sparse) data set, nchan x nchan x time x subjects x trials:
%-------------------------------------------------------------------------- 
if length(dims) == 5
    numtrials = size(mst_adj_in,5); 
    numsubj = size(mst_adj_in,4);
    numtimes = size(mst_adj_in,3);

    
end

%--------------------------------------------------------------------------        
% Full (non sparse) data set, nchan x nchan x time x subjects:
%-------------------------------------------------------------------------- 
if length(dims) == 4
    numsubj = size(mst_adj_in,4);
    numtimes = size(mst_adj_in,3);    
    numtrials = 1;
    
end

%--------------------------------------------------------------------------        
% No multiple subjects, nchan x nchan x time:
%-------------------------------------------------------------------------- 
if length(dims) < 4
    nsubj = 1;
    numtimes = size(mst_adj_in,3);
    ntrials = 1;

end

%--------------------------------------------------------------------------        
% No multiple timepoints, nchan x nchan:
%-------------------------------------------------------------------------- 
if length(dims) < 3 
    numsubj = 1;
    numtimes = 1; 
    numtrials = 1;
    
end

%--------------------------------------------------------------------------        
% Start loop across time points:
%-------------------------------------------------------------------------- 
for ntimes = 1:numtimes
    %ntimes = 1:size(mst_adj_in,3)

%--------------------------------------------------------------------------        
% Start loop across subjects:
%--------------------------------------------------------------------------             
    for nsubj = 1:numsubj
        %nsubj = 1:size(mst_adj_in,4)
        %mst_adj = mst_adj_in(:,:,ntimes,nsubj);
        %bin_adj = double(mst_adj ~= 0);

%--------------------------------------------------------------------------        
% Start loop across subjects:
%--------------------------------------------------------------------------             
        for ntrials = 1:numtrials
            %nsubj = 1:size(mst_adj_in,4)
            mst_adj = mst_adj_in(:,:,ntimes,nsubj,ntrials);
            bin_adj = double(mst_adj ~= 0);
        
%--------------------------------------------------------------------------        
% Strength (Nodal)
%--------------------------------------------------------------------------        
        GraphMeasuresOut.strength(:,ntimes,nsubj,ntrials) = SMILe_strengths_undirected(mst_adj);
        
%--------------------------------------------------------------------------        
% Degree (Nodal)
%--------------------------------------------------------------------------        
        GraphMeasuresOut.deg(:,ntimes,nsubj,ntrials) = SMILe_degrees_undirected(mst_adj);
        
%--------------------------------------------------------------------------
% Degree (Global)
%--------------------------------------------------------------------------        
        GraphMeasuresOut.avgdeg(ntimes,nsubj,ntrials) = SMILe_average_degree(mst_adj);
        
%--------------------------------------------------------------------------
% Global Hubs - TODO
%--------------------------------------------------------------------------                
        %[~, h] = SMILe_hub(mst_adj, deg, .05);
        %GraphMeasuresOut.hubs(ntimes,nsubj) = padcat(GraphMeasuresOut.hubs(ntimes,nsubj),h);
        %clear h
        
%--------------------------------------------------------------------------
% Eigenvalue Centrality (Nodal)
%--------------------------------------------------------------------------                
        GraphMeasuresOut.eig_centrality(:,ntimes,nsubj,ntrials) = SMILe_eigenvector_centrality_undirected(mst_adj);

%--------------------------------------------------------------------------
% Page Rank Centrality (Nodal)
%--------------------------------------------------------------------------                        
        r = SMILe_pagerank_centrality(mst_adj,0.85);
        GraphMeasuresOut.PageRank(:,ntimes,nsubj,ntrials) = r;
        
%--------------------------------------------------------------------------
% Betweenness Centrality (Nodal)
%--------------------------------------------------------------------------                
%         L = SMILe_weight_conversion(mst_adj, 'lengths');
%         BC = SMILe_betweenness_weighted(L);
%         BC = BC./[61*60];

        BC = SMILe_betweenness_binary(double(mst_adj ~= 0));
        BC = BC./[(size(mst_adj,1) - 1).*(size(mst_adj,1) - 2)];
        GraphMeasuresOut.BC(:,ntimes,nsubj,ntrials) = BC;
        clear L BC
        
%--------------------------------------------------------------------------
% Subgraph Centrality (Nodal)
%-------------------------------------------------------------------------- 
        SC = SMILe_subgraph_centrality(double(mst_adj ~= 0));
        GraphMeasuresOut.BC(:,ntimes,nsubj,ntrials) = SC;
        
%--------------------------------------------------------------------------                
% K-core Centrality
%--------------------------------------------------------------------------                
        [coreness,kn] = SMILe_kcoreness_centrality_binary_undirected(mst_adj);
        GraphMeasuresOut.kcore(:,ntimes,nsubj,ntrials) = coreness;

%--------------------------------------------------------------------------
% Edge Betweenness
%--------------------------------------------------------------------------
        [EBC, ~] = SMILe_edge_betweenness_binary(double(mst_adj ~= 0));
        GraphMeasuresOut.EBC(:,:,ntimes,nsubj,ntrials) = EBC;

%--------------------------------------------------------------------------
% Modularity (Global) and Community (Nodal)
%--------------------------------------------------------------------------                
        [Ci,Q] = SMILe_modularity_und(mst_adj, 3);
        GraphMeasuresOut.Community(:,ntimes,nsubj) = Ci;
        GraphMeasuresOut.Mod(ntimes,nsubj,ntrials) = Q;
        clear Ci Q

%--------------------------------------------------------------------------
% Modularity Degree Z-Score
%--------------------------------------------------------------------------                        
        Z = SMILe_module_degree_zscore(mst_adj,GraphMeasuresOut.Community(:,ntimes,nsubj),0);
        GraphMeasuresOut.ModDegZ(:,ntimes,nsubj,ntrials) = Z;
        
%--------------------------------------------------------------------------
% Leaf Fraction (Global)
%--------------------------------------------------------------------------
        GraphMeasuresOut.LF(ntimes,nsubj,ntrials) = SMILe_Leaf_Fraction(bin_adj);
        
%--------------------------------------------------------------------------
% Path Length (Global)
% Eccentricity (Nodal)
% Radius (Global)
% Diameter (Global)
%--------------------------------------------------------------------------                
        D = SMILe_distance_weighted(mst_adj);
        [lambda,~,ecc,radius,diameter] = SMILe_charpath(D, 0, 0);
        GraphMeasuresOut.PathLength(ntimes,nsubj,ntrials) = lambda;
        GraphMeasuresOut.Ecc(:,ntimes,nsubj,ntrials) = ecc;
        GraphMeasuresOut.Radius(ntimes,nsubj,ntrials) = radius;
        GraphMeasuresOut.Diameter(ntimes,nsubj,ntrials) = diameter;
        clear lambda ecc radius diameter
        
%--------------------------------------------------------------------------
% Efficiency (Nodal)
%--------------------------------------------------------------------------
        E = SMILe_efficiency(bin_adj,1);
        GraphMeasuresOut.Elocal(:,ntimes,nsubj) = E;
        clear E

%--------------------------------------------------------------------------
% Efficiency (Global)
%--------------------------------------------------------------------------
        E = SMILe_efficiency(bin_adj);
        GraphMeasuresOut.Eglobal(ntimes,nsubj,ntrials) = E;
        clear E
        
%--------------------------------------------------------------------------
% Transitivity (Global)
%--------------------------------------------------------------------------   
        T = SMILe_transitivity_weighted_undirected(mst_adj);
        GraphMeasuresOut.Radius(ntimes,nsubj,ntrials) = T;
        
%--------------------------------------------------------------------------
% Assortivity (Global)
%--------------------------------------------------------------------------                
        GraphMeasuresOut.Deg_Corr(ntimes,nsubj,ntrials) = SMILe_assortativity_binary(bin_adj,0);
        %clear bin_adj
        
%--------------------------------------------------------------------------
% Mean (Global)
%--------------------------------------------------------------------------                
        GraphMeasuresOut.mean(ntimes,nsubj,ntrials) = mean(mst_adj(mst_adj ~= 0));
        
%--------------------------------------------------------------------------
% Tree Hierarchy (Global)
%--------------------------------------------------------------------------                
        GraphMeasuresOut.Th(ntimes,nsubj,ntrials) = SMILe_TreeHierarchy(mst_adj);
        
%--------------------------------------------------------------------------
% Kappa (Global)
%--------------------------------------------------------------------------                
        %deg = SMILe_degrees_undirected(mst_adj);
        GraphMeasuresOut.K(ntimes,nsubj,ntrials) = SMILe_kappa(GraphMeasuresOut.deg(:,ntimes,nsubj));
        clear deg
        
%--------------------------------------------------------------------------
% Clustering Coefficient (Nodal)
%--------------------------------------------------------------------------
        if sum(mst_adj,1) == size(mst_adj,1)-1
            GraphMeasuresOut.CC(:,ntimes,nsubj,ntrials) = SMILe_clustering_coef_binary_undirected(mst_adj);

        else
            GraphMeasuresOut.CC(:,ntimes,nsubj,ntrials) = SMILe_clustering_coef_weighted_undirected(mst_adj);

        end
        
%--------------------------------------------------------------------------
% Clustering Coefficients (Global and Nodal)
%--------------------------------------------------------------------------        
        [C1,C2, C] = SMILE_clustCoeff(mst_adj);
        GraphMeasuresOut.GCC1(ntimes,nsubj,ntrials) = C1;
        GraphMeasuresOut.GCC2(ntimes,nsubj,ntrials) = C2;
        GraphMeasuresOut.NCC(:,ntimes,nsubj,ntrials) = C;
        
%--------------------------------------------------------------------------
% Density Measures (Global)
%--------------------------------------------------------------------------        
        [kden,N,NumConns] = SMILe_density_undirected(mst_adj);
        GraphMeasuresOut.kdensity(ntimes,nsubj,ntrials) = kden;
        GraphMeasuresOut.N(ntimes,nsubj,ntrials) = N;
        GraphMeasuresOut.NumConns(ntimes,nsubj,ntrials) = NumConns;
    
%--------------------------------------------------------------------------
% Small Worldness
%--------------------------------------------------------------------------
%        [S_ws, S_trans, S_ws_MC, S_trans_MC] = SMILe_SmallWorld(mst_adj);
        [S_ws, ~, ~, ~] = SMILe_SmallWorld(bin_adj);
        GraphMeasuresOut.SW(ntimes,nsubj,ntrials) = S_ws;
        
%--------------------------------------------------------------------------
% Reset MST
%--------------------------------------------------------------------------                
        mst_adj = [];
        bin_adj = [];

%--------------------------------------------------------------------------
% End for trials
%--------------------------------------------------------------------------                
        end
    
%--------------------------------------------------------------------------
% End for time subjects
%--------------------------------------------------------------------------            
    end
 
%--------------------------------------------------------------------------
% End for time bins
%--------------------------------------------------------------------------            
end

%--------------------------------------------------------------------------
% End of function
%--------------------------------------------------------------------------        
end