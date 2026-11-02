function [surr] = SMILe_Surrogates_RandTimePerm(datain, nsurr)
%--------------------------------------------------------------------------
% Random Permutation surrogates
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
% surr: The surrogate data sets in the same dimension as the input data
%-------------------------------------------------------------------------- 
%%
%-------------------------------------------------------------------------- 
% Get information about the size of the data set:
%-------------------------------------------------------------------------- 
origdims = size(datain);
ndims = numel(origdims);
if min(origdims) == 1
    ndims = 1;
    
end

if numel(origdims) == 3
    datain = squeeze(mean(datain,3));

end

dims = size(datain);
numchan = dims(1);

%-------------------------------------------------------------------------- 
% Allocate memory for output:
%-------------------------------------------------------------------------- 
surr = zeros([dims nsurr]);

%--------------------------------------------------------------------------
% Start progress bar:
%--------------------------------------------------------------------------
%wh = waitbar(0,'Starting...');

multiWaitbar( 'Channels', 'Color', [0.8 0.0 0.1] );
multiWaitbar( 'Surrogates', 'Color', [0.2 0.9 0.3] );

%--------------------------------------------------------------------------
% Loop across channels:
%--------------------------------------------------------------------------
for nchan = 1:numchan
    
%--------------------------------------------------------------------------
% Update waitbar:
%--------------------------------------------------------------------------
    multiWaitbar( 'Channels', nchan/numchan );
    
%-------------------------------------------------------------------------- 
% Reshape data to make life easy:
%--------------------------------------------------------------------------
    permdat = reshape(datain(nchan,:),1,[]);

%--------------------------------------------------------------------------
% Loop across surrogates:
%--------------------------------------------------------------------------
    for nsurr = 1:nsurr
        
%--------------------------------------------------------------------------
% Update waitbar:
%--------------------------------------------------------------------------        
        multiWaitbar( 'Surrogates', nsurr/nsurr );
%-------------------------------------------------------------------------- 
% Get a randomly permuted index of every element in permdat:
%--------------------------------------------------------------------------
        permOrder = randperm(numel(permdat));

%-------------------------------------------------------------------------- 
% Construct a permuted data set based on this index:
%--------------------------------------------------------------------------
        outsur(:,nsurr) = permdat(permOrder);

%-------------------------------------------------------------------------- 
% Recover original shape of the data:
%--------------------------------------------------------------------------
        %outsur = reshape(tmpsur,[dims]);

        
    end
        surr(nchan,:,:) = squeeze(outsur);
        
        
%-------------------------------------------------------------------------- 
% Gather surrogates:
%--------------------------------------------------------------------------
%     switch ndims
%         case 3
%             surr(:,:,:,nsurr) = outsur;
%         
%         case 2
%             surr(:,:,nsurr) = outsur;
%         
%         case 1
%             surr(:,nsurr) = outsur;
%         
%     end
    
end
 
%--------------------------------------------------------------------------
% Kill progress bar:
%--------------------------------------------------------------------------
%close(wh)
multiWaitbar('CloseAll');

%-------------------------------------------------------------------------- 
% End for function:
%--------------------------------------------------------------------------
end