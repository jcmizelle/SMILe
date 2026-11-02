function surr = SMILe_Shuffle(datain,numsurr)
%--------------------------------------------------------------------------
% Makes numsurr shuffled surrogates of a time series datain.
%--------------------------------------------------------------------------
% Inputs:
% datain: The original time series. Should be an N x 1 vector
% numsurr: The number of requested surrogates.
%--------------------------------------------------------------------------
% Output:
% surr: The shuffled time series in the form of time x numsurr. 

%%
%-------------------------------------------------------------------------- 
% Get information about the size of the data set:
%-------------------------------------------------------------------------- 
numpoints = length(datain);

%-------------------------------------------------------------------------- 
% Take care of missing input:
%-------------------------------------------------------------------------- 
if ~exist('numsurr', 'var') || isempty(numsurr)
   numsurr = 1;

end

%-------------------------------------------------------------------------- 
% Allocate memory:
%--------------------------------------------------------------------------
surr = zeros([numpoints numsurr]);

%-------------------------------------------------------------------------- 
% Loop across the number of requested surrogates:
%--------------------------------------------------------------------------
for i = 1:numsurr
    
%-------------------------------------------------------------------------- 
% Make a random vector with indices up to numpoints:
%--------------------------------------------------------------------------     
    tmpsurr = randperm(numpoints);
    
%-------------------------------------------------------------------------- 
% Shuffle the original datain with respect to the permutation vector:
%--------------------------------------------------------------------------         
    surr(:,i) = datain(tmpsurr);
    
end

%-------------------------------------------------------------------------- 
% End for function:
%-------------------------------------------------------------------------- 
end