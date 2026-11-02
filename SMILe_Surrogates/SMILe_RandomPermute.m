function [surr] = SMILe_RandomPermute(datain, nperms)
%--------------------------------------------------------------------------
% Random permutation surrogates of any N-dimensional data set.
%
% New data are created simply by random permutations of the original 
% series. The permutations guarantee the same amplitude distribution than 
% the original series, but destroy any linear correlation. This method is 
% associated to the null hypothesis of the data being uncorrelated noise 
% (possibly Gaussian and measured by a static nonlinear function.
%--------------------------------------------------------------------------
% Inputs:
% data: Typically a nchan x ntimes [x ntrials] matrix
% 
% nsurr: The number of desired surrogate data sets
%--------------------------------------------------------------------------
% Output:
% randout: The surrogate data sets in the same dimension as the input data
%-------------------------------------------------------------------------- 
%%
%--------------------------------------------------------------------------
% Get dimensions of the data:
%-------------------------------------------------------------------------- 
dims = size(datain);
ndims = numel(dims);
if min(dims) == 1
    ndims = 1;
    
end
% 
% if numel(dims) == 3
%     datain = squeeze(mean(datain,3));
% 
% end


%-------------------------------------------------------------------------- 
% Allocate memory for output:
%-------------------------------------------------------------------------- 
surr = zeros([dims nperms]);

%--------------------------------------------------------------------------
% Start progress bar:
%--------------------------------------------------------------------------
h = waitbar(0,'Starting...');

for nperm = 1:nperms
    waitbar(nperm/nperms, h, sprintf('Creating permutation %6.0f  of  %6.0f ', nperm, nperms))
%-------------------------------------------------------------------------- 
% Create matrix of randomly permuted integers to the total number of
% entries in the permutation matrix:
%-------------------------------------------------------------------------- 
    randidx = randperm(prod(dims));
    randmat = reshape(randidx,dims);

%-------------------------------------------------------------------------- 
% Reassign elements according to the randomization:
%-------------------------------------------------------------------------- 
    randout = datain(randmat);

%-------------------------------------------------------------------------- 
% Gather permutations:
%--------------------------------------------------------------------------

    switch ndims
        case 3
            surr(:,:,:,nperm) = randout;
        
        case 2
            surr(:,:,nperm) = randout;
        
        case 1
            surr(:,nperm) = randout;
        
    end
    
end

%--------------------------------------------------------------------------
% Kill progress bar:
%--------------------------------------------------------------------------
close(h)

%-------------------------------------------------------------------------- 
% End for function:
%--------------------------------------------------------------------------
end