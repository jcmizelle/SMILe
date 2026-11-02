function [cs, fftdataout, nave, nseg] = SMILe_CS_Fourier(datain, cfg)
%--------------------------------------------------------------------------
%
% This function calculates cross-spectra and for event-related and
% continuous measurements (with epoch-like event indicators)
%
% Input:
% data      - Typically a preprocessed data set (from the SMILe
%             preprocessing routine), but can be can be any data set that
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
% [cs, fftdataout, ~, ~] = SMILe_CS_Fourier(datain, cfg);
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
%$$$$$$ Add path(s) for Laplacian transform.
%$$$$$$
%$$$$$$ I always include this as a final strp and include the resulting
%$$$$$$ computations in the EEG structure.
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%addpath(genpath([eeglabdir filesep 'Laplace' filesep]));
%addpath(genpath('path/to/Laplace'));
%rehash


%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Compute the Laplacian transformation of the data and add to the EEG
%$$$$$$ structure, but not replacing the EEG.data structure.
%$$$$$$
%$$$$$$ Both measures are based on the Perrin method and yield similar
%$$$$$$ results. The CSD function also returns the spherical spline
%$$$$$$ surface potential (SP) interpolation.
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%[CSDout, SLAP] = laplacian_CSD(EEG);
%EEG.CSD = CSDout;
%EEG.SLAP = SLAP;
%data = EEG.CSD;

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
% Specify unset parameters.
%
% Check for input and set defaults. Assuming 1kHz sampling rate and an
% interest in frequencies less than 40Hz - if not true, set segleng to the
% actual sampling rate and desired frequency max.
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Sampling frequency - self explanatory. Assuming 1kHz. Real value.
%--------------------------------------------------------------------------
if ~exist('fs','var') || isempty(fs)
    fs = 1000;

end

%--------------------------------------------------------------------------
% Segment length. Defines frequency resolution. If segleng = fs, frequency
% resolution is 1Hz. Real value.
%--------------------------------------------------------------------------
if ~exist('segleng','var') || isempty(segleng)
    segleng = fs;

end

%--------------------------------------------------------------------------
% Amount of segment shift (overlapping or not). Real value.
%--------------------------------------------------------------------------
if ~exist('segshiftval','var') || isempty(segshiftval)
    segshiftval = 2;
    segshift = segleng/segshiftval;

end

%--------------------------------------------------------------------------
% Enforce a time delay per epoch/trial if not zero. Real value.
%--------------------------------------------------------------------------
if ~exist('segstart','var') || isempty(segstart)
    segstart = 0;

end

%--------------------------------------------------------------------------
% Type of window for later Fourier transform. Default is hanning as
% it works in practically every case. Real selection.
%--------------------------------------------------------------------------
if ~exist('wintype','var') || isempty(wintype)
    wintype = 'hanning';

end

%--------------------------------------------------------------------------
% Amount of segment shift (overlapping or not). Binary value (0 = don't
% zero pad, 1 = zero pad the data).
%--------------------------------------------------------------------------
if ~exist('zeropad','var') || isempty(zeropad)
    zeropad = 1;

end

%--------------------------------------------------------------------------
% Subtract average of the spectra. Binary value (0 = don't subtract the
% average, 1 = do subtract the average).
%--------------------------------------------------------------------------
if ~exist('subave','var') || isempty(subave)
    subave = 1;

end

%--------------------------------------------------------------------------
% Average across estimation segments. With segments, you kinda keep some
% idea of time resolution, but...
%
% Binary value (0 = don't average segments, 1 = average segments).
%--------------------------------------------------------------------------
if ~exist('segave','var') || isempty(segave)
    segave = 1;

end

%--------------------------------------------------------------------------
% Detrend data prior to Fourier. Binary value (0 = don't detrend,
% 1 = detrend).
%--------------------------------------------------------------------------
if ~exist('fftdetrend','var') || isempty(fftdetrend)
    fftdetrend = 1;

end

%--------------------------------------------------------------------------
% If selected, keep cross spectra per epoch. Adds an extra dimension to the
% output data. Binary value (0 = don't keep trials, 1 = keep trials).
%--------------------------------------------------------------------------
if ~exist('keeptrials','var') || isempty(keeptrials) || keeptrials == 0
    keeptrials = [];

end

%--------------------------------------------------------------------------
% If selected, keep cross spectra per epoch. Binary value (0 = don't keep
% trials, 1 = keep trials).
%--------------------------------------------------------------------------
if ~exist('maxfreq','var') || isempty(maxfreq)
    maxfreq = [];

end

%--------------------------------------------------------------------------
% The option to use a "band pass" for this process. Instead of pulling all
% frequencies to maxfreq, you could (e.g.) choose alpha by setting (in the
% cfg) band = [8 12]. This defaults to maxfreq. This also needs to be
% adjusted to accommodate bins and not real frequencies. This is handled
% later.
%
% Should be a vector.
%--------------------------------------------------------------------------
if ~exist('band','var') || isempty(band)
    band = [];
   
end


%--------------------------------------------------------------------------
% If data were not previously epoched (i.e., continuous), provide some
% value for repeated estimates. Guessing on appropriate number:
%--------------------------------------------------------------------------
if numel(origdims) == 3
    epleng = origdims(2);

else
    %epleng = floor(dims(1) / 50);
    epleng = dims(1);
   
end

%--------------------------------------------------------------------------
% Determine the number of eopchs and the total number of segments:
%--------------------------------------------------------------------------

segshiftval = 2;
segshift = segleng/segshiftval;
    
    
ndat = dims(1);
nchan = dims(2);
nep = floor(ndat/epleng);
nseg = floor((epleng-segleng)/segshift)+1;

%--------------------------------------------------------------------------
% Set windowing parameters:
%--------------------------------------------------------------------------
if strmatch('hanning', wintype)
    fftwindow = repmat(hanning(segleng),1,nchan);
   
elseif strmatch('hamming', wintype)
    fftwindow = repmat(hamming(segleng),1,nchan);
   
end

%--------------------------------------------------------------------------
% Zero pad if requested:
%--------------------------------------------------------------------------
if exist('zeropad','var') && ~isempty(zeropad) && zeropad ~= 0
    fftwindow = [fftwindow;zeros(2^nextpow2(segleng),nchan)];
   
end

%--------------------------------------------------------------------------
% Check for a band of interest. If not present, set a frequency range
% indicated as 0 Hz to maxfreq (these are in bins):
%--------------------------------------------------------------------------
if ~exist('band','var') | isempty(band) | band == 0
   
    if exist('maxfreq','var') && ~isempty(maxfreq)
        df = fs/segleng;
        maxfreqbin = round((maxfreq / df)) + 1;

    else
        df = fs/segleng;
        maxfreq = 40;
        maxfreqbin = round((maxfreq / df)) + 1;

    end
    disp('max freqency selected, and converted to apprpriate frequency bins')
    disp(maxfreqbin)

else
    disp('band selected -- converted to apprpriate frequency bins')
    disp(band)
    band = band+1;
   
end


%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Input and options / defaults prepared. Build output structures,
%$$$$$$ checking for options and defaults.
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

%--------------------------------------------------------------------------
% Prep output:
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% For 0 - max frequency range:
%--------------------------------------------------------------------------
if exist('maxfreqbin','var') && ~isempty(maxfreqbin)

%--------------------------------------------------------------------------
% If keeping trials:
%--------------------------------------------------------------------------    
    if exist('keeptrials', 'var') && ~isempty(keeptrials) && keeptrials > 0
       
        if segave == 0
            cs = zeros(nchan, nchan, maxfreqbin, nep, nseg);
            av = zeros(nchan, maxfreqbin, nep, nseg);

        else
            cs = zeros(nchan, nchan, maxfreqbin, nep);
            av = zeros(nchan, maxfreqbin, nep);

        end

%--------------------------------------------------------------------------
% If not keeping trials:
%--------------------------------------------------------------------------  
    else
        
        if segave == 0
            cs = zeros(nchan, nchan, maxfreqbin, nseg);
            av = zeros(nchan, maxfreqbin, nseg);

        else
            cs = zeros(nchan, nchan, maxfreqbin);
            av = zeros(nchan, maxfreqbin);

        end
       
    end
   
%--------------------------------------------------------------------------
% For a frequency band of interest:
%--------------------------------------------------------------------------    
elseif exist('band','var') && ~isempty(band)

%--------------------------------------------------------------------------
% If keeping trials:
%--------------------------------------------------------------------------    
    if exist('keeptrials', 'var') && ~isempty(keeptrials) && keeptrials > 0
       
        if segave == 0
            cs = zeros(nchan, nchan, nep, nseg);
            av = zeros(nchan, nep, nseg);

        else
            cs = zeros(nchan, nchan, nep);
            av = zeros(nchan, nep);

        end
       
%--------------------------------------------------------------------------
% If not keeping trials:
%--------------------------------------------------------------------------        
    else
        
        if segave == 0
            cs = zeros(nchan, nchan, nseg);
            av = zeros(nchan, nseg);

        else
            cs = zeros(nchan, nchan);
            av = zeros(nchan,1);

        end
       
    end
   
end

%--------------------------------------------------------------------------
% Build output for fft data:
%--------------------------------------------------------------------------
if exist('band','var') && ~isempty(band)
    
    if exist('keeptrials', 'var') && ~isempty(keeptrials) && keeptrials > 0
        
        if segave == 0
            fftdataout = zeros(nchan, nep, nseg);
            
        else
            fftdataout = zeros(nchan, nep);
        
        end
        
    else
        
        if segave == 0
            fftdataout = zeros(nchan, nseg);
            
        else
            fftdataout = zeros(nchan, 1);
        
        end
        
    end
   
else
    
    if exist('keeptrials', 'var') && ~isempty(keeptrials) && keeptrials > 0
        
        if segave == 0
            fftdataout = zeros(nchan, maxfreqbin, nep, nseg);
            
        else
            fftdataout = zeros(nchan, maxfreqbin, nep);

        end
        
    else
        
        if segave == 0
            fftdataout = zeros(nchan, maxfreqbin, nseg);
            
        else
            fftdataout = zeros(nchan, maxfreqbin);

        end
    end
       
end


%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Compute all Fourier estimates and compute cross spectra, paying
%$$$$$$ attention to all settings.
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

%--------------------------------------------------------------------------
% Keep track of the number of averages over which the cross spectrum is 
% estimated:
%--------------------------------------------------------------------------
nave = 0;

%--------------------------------------------------------------------------
% Start loop across epochs:
%--------------------------------------------------------------------------
for j = 1:nep

%--------------------------------------------------------------------------
% Waitbr for entertainment:
%--------------------------------------------------------------------------    
    waitbar(j/nep, h, sprintf('Working on trial: %6.0f ', j))

%--------------------------------------------------------------------------
% Now start:
%--------------------------------------------------------------------------    
    dataep = data((j-1)*epleng+1:j*epleng,:);
   
%--------------------------------------------------------------------------
% Start inner loop across segments within epochs:
%--------------------------------------------------------------------------
    for i = 1:nseg
        dataloc = dataep((i-1)*segshift+1:(i-1)*segshift+segleng,:);

        if fftdetrend == 1
            dataloc = detrend(dataloc);
           
        end

        if zeropad ~= 0
            dataloc = [dataloc;zeros(2^nextpow2(segleng),nchan)];
           
        end
       
%--------------------------------------------------------------------------
% Do the fft:
%--------------------------------------------------------------------------
        datalocfft = fft(dataloc.*fftwindow);
       
        if exist('band','var') && ~isempty(band)
            
            if exist('keeptrials', 'var') && ~isempty(keeptrials) && keeptrials > 0
                
                if segave == 0
                    fftdataout(:,j,i) = fftdataout(:,j,i) + mean(datalocfft(band(1):band(2),:),1)';
                    
                else
                    fftdataout(:,j) = fftdataout(:,j) + mean(datalocfft(band(1):band(2),:),1)';
           
                end
                
            else
                
                if segave == 0
                    fftdataout(:, i) = fftdataout(:, i) + mean(datalocfft(band(1):band(2),:),1)';
                    
                else
                    fftdataout(:) = fftdataout(:) + mean(datalocfft(band(1):band(2),:),1)';
           
                end
            
            end
                       
        else
            
            try 
                nfreq = maxfreqbin;
                
            catch
                nfreq = 41;
                
            end
            
            if exist('keeptrials', 'var') && ~isempty(keeptrials) && keeptrials > 0
                
                if segave == 0
                    fftdataout(:,1:nfreq,j,i) = fftdataout(:,1:nfreq,j,i) + datalocfft(1:nfreq,:)';
                    
                else
                    fftdataout(:,1:nfreq,j) = fftdataout(:,1:nfreq,j) + datalocfft(1:nfreq,:)';

                end
                
            else
                
                if segave == 0
                    fftdataout(:,1:nfreq,i) = fftdataout(:,1:nfreq,i) + datalocfft(1:nfreq,:)';
                    
                else
                    fftdataout(:,1:nfreq) = fftdataout(:,1:nfreq) + datalocfft(1:nfreq,:)';

                end
                
            end
        end

%--------------------------------------------------------------------------
% Compute cross spectra for a range of frequencies (0Hz - maxfreq):
%--------------------------------------------------------------------------
            if exist('maxfreqbin','var') && ~isempty(maxfreqbin)
               
                for f = 1:maxfreqbin
                    
%--------------------------------------------------------------------------
% Check keeptrials:
%--------------------------------------------------------------------------                  
                    if exist('keeptrials', 'var') && ~isempty(keeptrials) && keeptrials > 0
                       
                        if segave == 0
                            cs(:,:,f,j,i) = cs(:,:,f,j,i)+conj(datalocfft(f,:)'*datalocfft(f,:));
                            av(:,f,j,i) = av(:,f,j,i)+conj(datalocfft(f,:)');

                        else
                            cs(:,:,f,j) = cs(:,:,f,j)+conj(datalocfft(f,:)'*datalocfft(f,:));
                            av(:,f,j) = av(:,f,j)+conj(datalocfft(f,:)');

                        end
                        
%--------------------------------------------------------------------------
% Not keeping trials:
%--------------------------------------------------------------------------
                    else
       
                        if segave == 0
                            cs(:,:,f,i) = cs(:,:,f,i)+conj(datalocfft(f,:)'*datalocfft(f,:));
                            av(:,f,i) = av(:,f,i)+conj(datalocfft(f,:)');

                        else
                            cs(:,:,f) = cs(:,:,f)+conj(datalocfft(f,:)'*datalocfft(f,:));
                            av(:,f) = av(:,f)+conj(datalocfft(f,:)');

                        end

                    end
                   
                end
               
            end
               
%--------------------------------------------------------------------------
% Compute cross spectra for a specific band:
%--------------------------------------------------------------------------
            if exist('band','var') && ~isempty(band)
                
%--------------------------------------------------------------------------
% Check keeptrials:
%--------------------------------------------------------------------------                
                    if exist('keeptrials', 'var') && ~isempty(keeptrials) && keeptrials > 0
                       
                        if segave == 0
                            cs(:,:,j,i) = cs(:,:,j,i)+conj(mean(datalocfft(band(1):band(2),:),1)'*mean(datalocfft(band(1):band(2),:),1));
                            av(:,j,i) = av(:,j,i)+conj(mean(datalocfft(band(1):band(2),:),1)');

                        else
                            cs(:,:,j) = cs(:,:,j)+conj(mean(datalocfft(band(1):band(2),:),1)'*mean(datalocfft(band(1):band(2),:),1));
                            av(:,j) = av(:,j)+conj(mean(datalocfft(band(1):band(2),:),1)');

                        end
                        
%--------------------------------------------------------------------------
% Not keeping trials:
%--------------------------------------------------------------------------                      
                    else
       
                        if segave == 0
                            cs(:,:,i) = cs(:,:,i)+conj(mean(datalocfft(band(1):band(2),:),1)'*mean(datalocfft(band(1):band(2),:),1));
                            av(:,i) = av(:,i)+conj(mean(datalocfft(band(1):band(2),:),1)');

                        else
                            cs(:,:) = cs(:,:)+conj(mean(datalocfft(band(1):band(2),:),1)'*mean(datalocfft(band(1):band(2),:),1));
                            av(:) = av(:)+conj(mean(datalocfft(band(1):band(2),:),1)');
                           
                        end

                    end
                   
            end
   
%--------------------------------------------------------------------------
% End loop across segments:
%--------------------------------------------------------------------------
    end
   
%--------------------------------------------------------------------------
% Update counter:
%--------------------------------------------------------------------------
nave = nave+1;

%--------------------------------------------------------------------------
% End loop across epochs:
%--------------------------------------------------------------------------
end

%--------------------------------------------------------------------------
% Kill progress bar:
%--------------------------------------------------------------------------
close(h)


%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Fourier estimation and cross spectra computed. Normalize them as a
%$$$$$$ function of the number of averages and / or segments depending on
%$$$$$$ settings.
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
if segave == 0
    cs = cs/nave;
    av = av/nave;
    fftdataout = fftdataout/nave;

else
    nave = nave*nseg;  
    cs = cs/nave;
    av = av/nave;
    fftdataout = fftdataout/nave;

end


%%
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%$$
%$$$$$$ Check for subtracting the average spectrum. Do it if found.
%$$
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
if subave == 1

%--------------------------------------------------------------------------
% Check for a range of frequency range:
%--------------------------------------------------------------------------
    if exist('maxfreqbin','var') && ~isempty(maxfreqbin)

%--------------------------------------------------------------------------
% Check for keeptrials:
%--------------------------------------------------------------------------      
        if exist('keeptrials', 'var') && ~isempty(keeptrials) && keeptrials > 0
       
             for f = 1:maxfreqbin
             
                if segave == 0
             
                    for i = 1:nseg
                       
                        for j = 1:nep
                            cs(:,:,f,j,i) = cs(:,:,f,j,i)-av(:,f,j,i)*av(:,f,j,i)';
                           
                        end
                   
                    end
                   
                else
                   
                    for j = 1:nep
                        cs(:,:,f,j) = cs(:,:,f,j)-av(:,f,j)*av(:,f,j)';
                   
                    end
                   
                end

            end
            
%--------------------------------------------------------------------------
% Not keeping trials:
%--------------------------------------------------------------------------
        else
           
            for f = 1:maxfreqbin
               
                if segave == 0
                   
                    for i = 1:nseg
                    cs(:,:,f,i) = cs (:,:,f,i)-av(:,f,i)*av(:,f,i)';

                    end

                else
                    cs(:,:,f) = cs(:,:,f)-av(:,f)*av(:,f)';

                end
               
            end
           
        end
       
    end
   
%--------------------------------------------------------------------------
% Specific band:
%--------------------------------------------------------------------------
elseif exist('band','var') && ~isempty(band)

%--------------------------------------------------------------------------
% Check for keeptrials:
%--------------------------------------------------------------------------
    if exist('keeptrials', 'var') && ~isempty(keeptrials) && keeptrials > 0

        if segave == 0
           
            for i = 1:nseg
                cs(:,:,j,i) = cs (:,:,j,i)-av(:,j,i)*av(:,j,i)';

            end

        else
            cs(:,:,j) = cs(:,:,j)-av(:,j)*av(:,j)';

        end
        
%--------------------------------------------------------------------------
% Not keeping trials:
%--------------------------------------------------------------------------          
    else
        if segave == 0
            cs(:,:,i) = cs (:,:,i)-av(:,i)*av(:,i)';

        else
            cs(:,:,j) = cs(:,:,j)-av(:,j)*av(:,j)';
           
        end
           
    end
       
end

%--------------------------------------------------------------------------
% End for parent function:
%--------------------------------------------------------------------------
end