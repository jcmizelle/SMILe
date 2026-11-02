function [cs, fftdataout, nave, nseg] = SMILe_CS_Fourier_Surrogate(datain, cfg)
%--------------------------------------------------------------------------
%
% This function calculates cross-spectra and for event-related and
% continuous measurements (with epoch-like event indicators)
%
% Input:
% data      - Typically a proprocessed data set (from the SMILe
%             proprocessing routine), but can be can be any data set that
%             is either in nchan x time x epochs or nchan x time format.
%             Reshaping is automatic.
%
% cfg       - A collection of configuration parameters (identified and 
%             defined below) that specify much of the internal mechanics 
%             of this file. If you are unsure about your data, leave 
%             everything in cfg set to [], and let the internal defaults 
%             set. An example to strictly follow is below.
%
% Output:
% nave:     - The number of averages over which the cross spectrum was
%             estimated
%
% cs:       - The shape of the output directly depends on your input data
%             and the various configurations that are selected. For 
%             example, if the "band" option is selected in addition to the 
%             "segave" option, then cs = nchan x nchan. If the same is 
%             true, but you do not average across segments, then 
%             cs = nchan x nchan x nseg
%
%           - If "maxfreq" is set, then cs = nchan x nchan x numfreq, and
%             further if you do not average across segments, then
%             cs = nchan x nchan x numfreq x nseg
%
%           - It can get complicated. In any case:
%
%             cs(nchan,nchan) will always be true. Additional tensors 
%             include frequencies (0 - maxfreq if not band selected), and 
%             epochs (if selected as an output for later connectivity 
%             calculations)
%
%             Possible output dimensions for cs:
%
%             cs = nchan x nchan (band)
%             cs = nchan x nchan (band) x segment
%             cs = nchan x nchan (band) x epoch
%             cs = nchan x nchan (band) x epoch x segment
%             cs = nchan x nchan x freqs x segment
%             cs = nchan x nchan x freqs x epoch
%             cs = nchan x nchan x freqs x epoch x segment
%
%
%--------------------------------------------------------------------------
% Example usage:
% clear cs nave cfg
%
% cfg = [];
% cfg.segleng = [];
% cfg.segshiftval = [];
% cfg.segstart = [];
% cfg.fs = [];
% cfg.wintype = [];
% cfg.zeropad = [];
% cfg.subave = [1];
% cfg.segave = [1];
% cfg.fftdetrend = [];
% cfg.keeptrials = [];
% cfg.maxfreq = [];
% cfg.band = [8 12];
% 
% [cs, fftdataout] = SMILe_CS_Fourier(datain, cfg);
%
%--------------------------------------------------------------------------


%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Get started
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

%--------------------------------------------------------------------------
% Start progress bar:
%--------------------------------------------------------------------------
h = waitbar(0,'Starting...');

%--------------------------------------------------------------------------
% Check for structure:
%--------------------------------------------------------------------------
if isstruct(datain)
    data = datain.data;
    
else 
    data = datain;
    
end

%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Bring data in and start checking dimensions, check parameters and
%$$$$$$ set defaults if options are not present.
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

%--------------------------------------------------------------------------
% Check size of data. Assume either nchan x time x epochs (EEGLAB format 
% that is the default output for our lab) or pre-formatted data in 
% time x nchan:
%--------------------------------------------------------------------------
origdims = size(data);

%--------------------------------------------------------------------------
% Convert to 2-D if in 3-D:
%--------------------------------------------------------------------------
if numel(origdims) > 2
    data = reshape(data,[origdims(1),origdims(2)*origdims(3)]);

end

%--------------------------------------------------------------------------
% Permute to time x nchan if not in that format and reset dims:
%--------------------------------------------------------------------------
dims = size(data);

if dims(2) > dims(1)
    data = permute(data, [2 1]);

end

dims = size(data);

%--------------------------------------------------------------------------
% Unpack configuration options, and assign parameter fields and values to
% variables in the local workspace:
%--------------------------------------------------------------------------
fn = fieldnames(cfg);

for nfield = 1:numel(fn)
    fni = char(fn(nfield));
    field = cfg.(fni);

    feval(@()assignin('caller', fni, field));
    %feval(@()assignin('base', fni, field));

end

%--------------------------------------------------------------------------
% Create surrogate data from real data:
%--------------------------------------------------------------------------

for nperm = 1:1000;
randtimes = randperm(size(data,1));
randidx = randperm(size(data,2));

surrogate_data(1:dims(1),1:dims(2)) = data(randtimes,randidx);

clear data;
data = surrogate_data;

[cs, ~, nave, nseg] = SMILe_CS_Fourier(EEG, cfg);
csout(:,:,nperm) = cs;

%--------------------------------------------------------------------------
% End for number of permutations:
%--------------------------------------------------------------------------
end

%--------------------------------------------------------------------------
% End for parent function:
%--------------------------------------------------------------------------
end