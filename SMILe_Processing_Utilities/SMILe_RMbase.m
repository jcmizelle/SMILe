function [dataout, datamean] = SMILe_RMbase(data, frames, basevector)
%% 
% Usage: 
% To remove whole-data channel means:
% [dataout] = rmbase(data); 
% 
% To remove mean of basevector from each channel and epoch:
% [dataout datamean] = rmbase(data, frames, basevector);
%            
% Input:
% data       - data matrix (chans, frames*epochs) or (chans, frames, epochs);
% frames     - data points per epoch {[]|0|default->data length}
% basevector - vector of baseline frames per epoch
%                 Ex 1:128 {[]|0|default->whole epoch}
%
%%
%--------------------------------------------------------------------------
% Start parsing input:
%--------------------------------------------------------------------------
if nargin < 3
    basevector = 0;
    
end

if isempty(basevector)
    basevector = 0;
    
end

if length(basevector) == 1 && basevector(1) ~= 0
   fprintf('SMILe_RMbase(): basevector should be a vector of frame indices.\n')
   
   return

end

if sum(basevector < 0) ~= 0
   fprintf('SMILe_RMbase(): basevector should be 0 or a vector of positive frame indices.\n');
   
   return

end

if nargin < 2
    frames = 0;
    
end

if isempty(frames)
    frames = 0;
    
end

if nargin < 1
    fprintf('SMILe_RMbase(): needs at least one argument.\n\n');
    
    return
    
end

if isempty(data)
    fprintf('SMILe_RMbase(): input data is empty.\n\n');
    
    return
    
end

%%
%--------------------------------------------------------------------------
% Handle data dimensions and shape:
%--------------------------------------------------------------------------
origdims = size(data);

if origdims(data) == 3
    data = reshape(data, size(data,1), size(data,2)*size(data,3));
    do_reshape = 1;
    
end	

[nchans, ntimes]= size(data);
if frames == 0
    frames = ntimes;

end  

    
epochs = fix(ntimes/frames);

if length(basevector)>ntimes
    fprintf('SMILe_RMbase(): length(basevector) > frames per epoch.\n\n');
    
    return
    
end

%--------------------------------------------------------------------------
% Set output data:
%--------------------------------------------------------------------------
datamean = zeros(nchans, epochs);
dataout = data;

%--------------------------------------------------------------------------
% Handle calculations:
%--------------------------------------------------------------------------
for e = 1:epochs
    
    for c = 1:nchans
        
        if basevector(1) ~= 0
            rmeans = nan_mean(double(data(c,(e-1)*frames+basevector)'));
            
        else
            rmeans = nan_mean(double(data(c,(e-1)*frames+1:e*frames)'));

        end
        
        datamean(c,e) = rmeans;
        dataout(c,(e-1)*frames+1:e*frames) = data(c,(e-1)*frames+1:e*frames) - rmeans;
    
    end
    
end

%--------------------------------------------------------------------------
% Reshape if required:
%--------------------------------------------------------------------------
if exist('do_reshape', 'var') && do_reshape == 1
dataout = reshape(dataout, origdims);
   
end

%--------------------------------------------------------------------------
% End of function
%--------------------------------------------------------------------------
end
