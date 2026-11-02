function adj = SMILe_Simulate_Conn_Trials(cs, nsim, pctdat, type)
%%
%--------------------------------------------------------------------------
% 
% This function returns an nchan x nchan x reps matrix of simulated
% adjacency matrices that can be used for further processing (e.g., graph
% theory, pattern classification, etc.).
%
% This is effectively a Monte Carlo style of simulation, in which simulated 
% trials are generated through random sampling of "truth" data to create 
% resampled data trials (or participants, etc.).
%
%--------------------------------------------------------------------------
%
% Mandatory input:
% - cs: complex valued cross spectrum of dimensions nchan x nchan x trials.
%   This should have come from one of the SMILe codes that calculates the 
%   cross spectrum (Fourier or wavelet)
%
% Optional input:
% - nsim: how many simulated trials are required?
%
% - pctdat: what percentage of the original trials should be used for each
%   (random) permutated simulation?
% 
% - type: which connectivity calculation should be used?
%
%--------------------------------------------------------------------------

%%
%--------------------------------------------------------------------------
% Set the number of simulations:
%--------------------------------------------------------------------------
if ~exist('nsim', 'var') || isempty(nsim)
    nsim = 50;
    
else
    nsim = nsim;
    
end

%--------------------------------------------------------------------------
% Set the percentage of data to use per simulation:
%--------------------------------------------------------------------------
if ~exist('pctdat', 'var') || isempty(pctdat)
    pctdat = .6;
    
elseif exist('pctdat', 'var') && pctdat > 1
    pctdat = pctdat / 100;
    
else
    pctdat = pctdat;
    
end

%--------------------------------------------------------------------------
% Set the default connectivity type:
%--------------------------------------------------------------------------
if ~exist('type', 'var') || isempty(type)
    type = 'wPLI';
    
else
    type = type;
    
end

%--------------------------------------------------------------------------
% Preallocate memory:
%--------------------------------------------------------------------------
adj = zeros(size(cs,1), size(cs,2), nsim);

%--------------------------------------------------------------------------
% Loop across desired simulations:
%--------------------------------------------------------------------------
for iter = 1:nsim
    
%--------------------------------------------------------------------------
% Generate random sequence to draw truth trials from:
%--------------------------------------------------------------------------    
    perms = randperm(size(cs,3));
    
%--------------------------------------------------------------------------
% Do the desired calculation - switch/case:
%--------------------------------------------------------------------------    
    switch type
        case 'PLI'
            adj(:,:,iter) = SMILe_PLI(squeeze(cs(:,:,perms(1:ceil(size(perms,2)*pctdat)))));
            
        case 'wPLI'
            adj(:,:,iter) = SMILe_wPLI(squeeze(cs(:,:,perms(1:ceil(size(perms,2)*pctdat)))));
    
        case 'dsPLI'
            adj(:,:,iter) = SMILe_dsPLI(squeeze(cs(:,:,perms(1:ceil(size(perms,2)*pctdat)))));
            
        case 'dswPLI'
            adj(:,:,iter) = SMILe_dswPLI(squeeze(cs(:,:,perms(1:ceil(size(perms,2)*pctdat)))));
        
        case 'PLV'
            adj(:,:,iter) = SMILe_PLV(squeeze(cs(:,:,perms(1:ceil(size(perms,2)*pctdat)))));                        
            
        case 'iPLV'
            adj(:,:,iter) = SMILe_iPLV(squeeze(cs(:,:,perms(1:ceil(size(perms,2)*pctdat)))));
            
        case 'ciPLV'
            adj(:,:,iter) = SMILe_ciPLV(squeeze(cs(:,:,perms(1:ceil(size(perms,2)*pctdat)))));
%--------------------------------------------------------------------------
% End for switch/case:
%--------------------------------------------------------------------------            
    end
    
%--------------------------------------------------------------------------
% End for the number of simulations:
%--------------------------------------------------------------------------    
end

%--------------------------------------------------------------------------
% End for functino
%--------------------------------------------------------------------------
end