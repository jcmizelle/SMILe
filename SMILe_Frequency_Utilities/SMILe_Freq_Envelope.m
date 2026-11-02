function [envelope_bounds] = SMILe_Freq_Envelope(fdata, npts, type)
%%
%--------------------------------------------------------------------------
% 
% This function returns the envelope of a signal. Here, we take advantage
% of the Matlab built-in function "envelope.m" and have access to three
% different envelope types: Analytic, RMS and Peak(s). This should be
% specified in the input.
%
%--------------------------------------------------------------------------
%
% Mandatory input:
% - fdata: this sould have come from the SMILe filter functions and
%   represents the band-pass frequency resolved signal. Can be a vector, a
%   2-D matrix of samples x channels or a or 3-D matrix of 
%   samples x channels x trials (or participants, etc.).
%
% Optional input:
% - npts: specifies the length of the filter (Hilbert) in the 'Analytic' 
%   option, number of points used in the 'RMS' option, or the separation of 
%   peaks in the 'Peaks' option. Default is 30.
%
% - type: method of envelope calculation. Can be either 'analytic', 'rms'
%   or 'peak'. Default is 'peak'
%
% Display of function
% fdatain = fdata(1,:,1);
% [envelope_bounds] = SMILe_Freq_Envelope(fdatain, 30, 'peak');
% plot(envelope_bounds.upper)
% hold on
% plot(envelope_bounds.lower,'k')
% plot(fdata(1,:,1),'r')
%
% figure;
% x = linspace(1,750,750);
% y1 = envelope_bounds.upper(17,:,5);
% y2 = envelope_bounds.lower(17,:,5);
% SMILe_Shade(x,y1,x,y2,'LineWidth',2, 'FillType',[1 2;2 1], 'FillColor',[.5 .5 .5]);
% SMILe_Shade(x,y1,'LineWidth',2,'FillType',[1 0; 0 1]);
% SMILe_Shade(x,y1,x,y2+2,'LineWidth',2,'FillType',[1 2; 1 0]);
% 
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Set the default filter length:
%--------------------------------------------------------------------------
if ~exist('flength', 'var') || isempty(npts)
    npts = 10;
    
else
    npts = npts;
    
end

%--------------------------------------------------------------------------
% Set the default envelope type:
%--------------------------------------------------------------------------
if ~exist('type', 'var') || isempty(type)
    type = 'peak';
    
else
    type = type;
    
end
%--------------------------------------------------------------------------
% Get dimensions of the data, use switch/case for dimensions:
%--------------------------------------------------------------------------
dims = size(fdata);
numdims = numel(dims);
if min(dims) == 1
    numdims = 1;
    
end

switch numdims
    
    case 3
        
%--------------------------------------------------------------------------
% Loop across trials:
%--------------------------------------------------------------------------
        for ntrial = 1:dims(3)
   
%--------------------------------------------------------------------------
% Loop across channels and calculate the envelope:
%--------------------------------------------------------------------------
            for nchan = 1:dims(1)
                [yup ylow] = envelope(squeeze(fdata(nchan,:,ntrial)), npts, type);
                envUp(nchan, :, ntrial) = yup;
                envLow(nchan, :, ntrial) = ylow;

%--------------------------------------------------------------------------
% End channel loop:
%--------------------------------------------------------------------------
            end
            
%--------------------------------------------------------------------------
% End trial loop:
%--------------------------------------------------------------------------
        end

    case 2   
        
%--------------------------------------------------------------------------
% Loop across channels and calculate the envelope:
%--------------------------------------------------------------------------
            for nchan = 1:dims(1)
                [yup ylow] = envelope(squeeze(fdata(nchan,:)), npts, type);
                envUp(nchan, :) = yup;
                envLow(nchan, :) = ylow;

%--------------------------------------------------------------------------
% End channel loop:
%--------------------------------------------------------------------------
            end
        
    case 1
        
%--------------------------------------------------------------------------
% Single channel/trial/signal, calculate the envelope:
%--------------------------------------------------------------------------
            [yup ylow] = envelope(fdata, npts, type);
            envUp = yup;
            envLow = ylow;
   
%--------------------------------------------------------------------------
% End for switch/case:
%--------------------------------------------------------------------------        
end

%--------------------------------------------------------------------------
% Collect output data:
%--------------------------------------------------------------------------
envelope_bounds.upper = envUp;
envelope_bounds.lower = envLow;

%--------------------------------------------------------------------------
% End for function:
%--------------------------------------------------------------------------
end