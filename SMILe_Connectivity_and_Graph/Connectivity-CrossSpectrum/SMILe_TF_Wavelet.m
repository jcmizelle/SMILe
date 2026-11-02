function [crss, powerspec, connmat, itpc, params] = SMILe_TF_Wavelet(datain, cfg)
%--------------------------------------------------------------------------
%
% This function calculates the multichannel time-frequency properties for a
% given input data set.
%
% Requires the function demean, in SMILe Utilites on the Z drive
%
% Input:
% datain: Preprocesed data set from our SMILe routine. Needs to be epoched.
%
% cfg: Processing parameters as outlined here
%       cfg = [];
%       cfg.savetime =      -200:10:1250; in ms, time 0 is stimulus onset
%       cfg.basetime =      [-1500 -750]; in ms, time 0 is stimulus onset
%       cfg.range_cycles =  [ 4 10 ]; leave default unless you have a
%                           reason to change
%       cfg.min_freq =  1;  minimum frequency of interest
%       cfg.max_freq = 45;  maximum frequency of interest
%       cfg.num_frex = 45;  number of frequency steps between min and max
%       cfg.non_phase = 0;  do / do not compute non phase locked values
%
% Output:
%--------------------------------------------------------------------------
% The actual cross spectrum to be used in other connectivity analyses:
%--------------------------------------------------------------------------
% crss.crssspctrm
% crss.blcrssspctrm
%
%--------------------------------------------------------------------------
% Output related to spectral power variables:
%--------------------------------------------------------------------------
% powerspec.tf_pow_trials = tf_pow_trials;
% powerspec.tf_pow_mean = tf_pow_mean;
% powerspec.tf_pow_baseline = tf_pow_baseline;
% powerspec.tf_pow_db = tf_pow_db;
% powerspec.tf_pow_pc = tf_pow_pc;
%
%--------------------------------------------------------------------------
% Output related to connectivity:
%--------------------------------------------------------------------------
% connmat.pli = pli;
% connmat.ispc = ispc;
% connmat.blpli = blpli;
% connmat.blispc = blispc;
%
%--------------------------------------------------------------------------
% Inter trial phase clustering:
%--------------------------------------------------------------------------
% itpc
%
%--------------------------------------------------------------------------
%
%--------------------------------------------------------------------------
% Information related to the parameters:
%--------------------------------------------------------------------------
% params.datatimes = datatimes;
% params.frex = frex;
% params.savetime = savetime;
% params.basetime = basetime;
% params.range_cycles = range_cycles;
% params.max_freq = max_freq;
% params.min_freq = min_freq;
%
%--------------------------------------------------------------------------
% Check for structure:
%--------------------------------------------------------------------------
if isstruct(datain)
    data = datain.data;
   
    fn = fieldnames(datain);

    for nfield = 1:numel(fn)
        fni = char(fn(nfield));
        field = datain.(fni);

        feval(@()assignin('caller', fni, field));
        %feval(@()assignin('base', fni, field));

    end
   
else
    data = datain;
   
end

fn = fieldnames(cfg);

for nfield = 1:numel(fn)
    fni = char(fn(nfield));
    field = cfg.(fni);

    feval(@()assignin('caller', fni, field));
    %feval(@()assignin('base', fni, field));

end

%--------------------------------------------------------------------------
% Get size of the data:
%--------------------------------------------------------------------------
dims = size(data);

if ~exist('trials', 'var') || isempty(trials)
    trials = dims(3);
   
end

if isfield(datain, 'times')%~exist('times', 'var') || isempty(times)
    datatimes = datain.times;
   
else
    datatimes = times;
   
end

%--------------------------------------------------------------------------
% Settings for times-of-interest and also for the baseline normalization
% (for power):
%--------------------------------------------------------------------------
if ~exist('savetime', 'var') || isempty(savetime)
    savetime = -200:10:1250;
   
end

if ~exist('basetime', 'var') || isempty(basetime)
    basetime = [datatimes(1) datatimes(1) + 250];
   
end

%--------------------------------------------------------------------------
% Get index of selected times:
%--------------------------------------------------------------------------
savetimeidx = zeros(size(savetime));

for ti = 1:length(savetimeidx)
    [~, savetimeidx(ti)] = min(abs(datatimes-savetime(ti)));
   
end

%--------------------------------------------------------------------------
% Get index of baseline times:
%--------------------------------------------------------------------------
if sum(abs(basetime))
[~, basetimeidx(1)] = min(abs(datatimes-basetime(1)));
    [~, basetimeidx(2)] = min(abs(datatimes-basetime(2)));
       
end

%--------------------------------------------------------------------------
% Other wavelet parameters. Temporal and frequency precision is governed by
% the number of cycles used in the wavelet creation. Here, we will use a
% range of cycles to hopefully get both good temporal precision and good
% frequency precision.
%
% Check input and set defaults if not found:
%--------------------------------------------------------------------------
if ~exist('range_cycles', 'var') || isempty(range_cycles)
    range_cycles = [ 4 10 ];
   
end

if ~exist('min_freq', 'var') || isempty(min_freq)
    min_freq =  1;
   
end

if ~exist('max_freq', 'var') || isempty(max_freq)
    max_freq = 45;
   
end

if ~exist('num_frex', 'var') || isempty(num_frex)
    %num_frex = 44;
    num_frex = (max_freq - min_freq) + 1;
   
end

if ~exist('non_phase', 'var') || isempty(non_phase)
    non_phase = 0;
   
end

%--------------------------------------------------------------------------
% Build a frequency vector:
%--------------------------------------------------------------------------
frex = linspace(min_freq, max_freq, num_frex);

%--------------------------------------------------------------------------
% Create a vector of wavelet cycle widths:
%--------------------------------------------------------------------------
wavecycles = linspace(range_cycles(1), range_cycles(end), num_frex) ./ (2*pi*frex);

%--------------------------------------------------------------------------
% Set up wavelet time length:
%--------------------------------------------------------------------------
wavtime = -3:1/srate:3;

%--------------------------------------------------------------------------
% Only use half of the spectrum (and wavelet):
%--------------------------------------------------------------------------
half_wave = (length(wavtime)-1)/2;

%--------------------------------------------------------------------------
% Length of the wavelet:
%--------------------------------------------------------------------------
nWave = length(wavtime);

%--------------------------------------------------------------------------
% Number of Fourier data points. Below, we will concatenate all trials, so we need
% the total number of points (all times across all trials:
%--------------------------------------------------------------------------
nData = pnts * trials;

%--------------------------------------------------------------------------
% Total length of the FFT:
%--------------------------------------------------------------------------
nConv = nWave + nData - 1;

%--------------------------------------------------------------------------
% Preallocate for time-frequency data. This is a situation where speed is
% dramatically increased:
%--------------------------------------------------------------------------
full_analytic_signal = zeros(dims(1), size(frex, 2), pnts, trials);
full_pow = zeros(dims(1), size(frex, 2), pnts, trials);
bins_analytic_signal = zeros(dims(1), size(frex, 2), size(savetimeidx,2), trials);
bins_pow = zeros(dims(1), size(frex, 2), size(savetimeidx,2), trials);
baseline_analytic_signal = zeros(dims(1), size(frex, 2), trials);
baselinepow = zeros(dims(1), size(frex, 2), trials);

%--------------------------------------------------------------------------
% Fancy waitbars for entertainment and to monitor progress:
%--------------------------------------------------------------------------
multiWaitbar( 'Channels', 'Color', [0.8 0.0 0.1] );
multiWaitbar( 'Frequencies', 'Color', [0.2 0.9 0.3] );
multiWaitbar( 'Computing Spectral Properties - Channel j', 'Color', [1.0 0.4 0.0] );
multiWaitbar( 'Computing Spectral Properties - Channel i', 'Color', [0.1 0.5 0.8] );

%--------------------------------------------------------------------------
% Start loop across channels:
%--------------------------------------------------------------------------  
for nchan = 1:dims(1)
%--------------------------------------------------------------------------
% Update waitbar:
%--------------------------------------------------------------------------
    multiWaitbar( 'Channels', nchan/dims(1) );

%--------------------------------------------------------------------------
% Optional ERP subraction (demean) to compute non-phase locked power:
%--------------------------------------------------------------------------
    if exist('non_phase', 'var') && ~isempty(non_phase) && non_phase > 0
        dimsize = size(data,3);
        dimrep = ones(1,length(dims));
        dimrep(3) = dimsize;
        data = data - repmat(mean(data, 3), dimrep);

    end

%--------------------------------------------------------------------------
% Concatenate trials:
%--------------------------------------------------------------------------
    alldata = reshape(data(nchan,:,:), 1, []);
    alldata = detrend(alldata);

%--------------------------------------------------------------------------
% Compute FFT of the data. This is independent of trials and frequency
% (spectrum of concatenated data will maintain the properties of the
% individual trials along its length):
%--------------------------------------------------------------------------
    dataX = fft( alldata ,nConv );

%--------------------------------------------------------------------------
% Start loop across frequencies:
%--------------------------------------------------------------------------  
    for fi = 1:length(frex)

%--------------------------------------------------------------------------
% Update waitbar:
%--------------------------------------------------------------------------
        multiWaitbar( 'Frequencies', fi/length(frex) );

%--------------------------------------------------------------------------
% Create wavelet and get its FFT. The wavelet doesn't change on each trial,
% so we only need to calculate this once:
%--------------------------------------------------------------------------
        wavelet  = exp(2*1i*pi*frex(fi).*wavtime) .* exp(-wavtime.^2./(2*wavecycles(fi)^2));
        waveletX = fft(wavelet,nConv);
        waveletX = waveletX ./ max(waveletX);

%--------------------------------------------------------------------------
% Now, run convolution in one step to get our analytic signal:
%--------------------------------------------------------------------------
        as = ifft(waveletX .* dataX);
        as = as(half_wave+1:end-half_wave);

%--------------------------------------------------------------------------
% Reshape back to time x trials:
%--------------------------------------------------------------------------
        as = reshape( as, pnts, trials );
     
%--------------------------------------------------------------------------
% Collect full analytic signal:
%--------------------------------------------------------------------------
        full_analytic_signal(nchan, fi, :, :) = as;

%--------------------------------------------------------------------------
% Collect the full power spectrum:
%--------------------------------------------------------------------------
        full_pow(nchan, fi, :, :) = abs(as).^2;

%--------------------------------------------------------------------------
% Downsample - collect the analytic signal time points of interest:
%--------------------------------------------------------------------------
        bins_analytic_signal(nchan, fi, :, :) = as(savetimeidx, :);
       
%--------------------------------------------------------------------------
% Downsample - collect the power spectrum at time points of interest:
%--------------------------------------------------------------------------
        bins_pow(nchan, fi, :, :) = abs(as(savetimeidx, :).^2);

%--------------------------------------------------------------------------
% Collect analytic signal baseline points:
%--------------------------------------------------------------------------
        baseline_analytic_signal(nchan, fi, :) = mean(as(basetimeidx(1):basetimeidx(2),:), 1);

%--------------------------------------------------------------------------
% Collect power spectrum baseline points:
%--------------------------------------------------------------------------
        baselinepow(nchan, fi, :) = mean(abs(as(basetimeidx(1):basetimeidx(2),:)).^2, 1);

%--------------------------------------------------------------------------
% End for frequencies:
%--------------------------------------------------------------------------      
    end

%--------------------------------------------------------------------------
% End for channels:
%--------------------------------------------------------------------------      
end

%--------------------------------------------------------------------------
% Power spectrum, all time points and trials:
%--------------------------------------------------------------------------
tf_pow_trials = abs(full_analytic_signal.^2);

%--------------------------------------------------------------------------
% Power spectrum, all time points - averaged across trials:
%--------------------------------------------------------------------------
tf_pow_mean = mean(abs(bins_analytic_signal).^2,4);
% tf_pow_bl = 10*log10( tf_pow ./ baselinedata );
%tf_pow_bl = 10*(log10(tf_pow) - log10(baselinedata));

%--------------------------------------------------------------------------
% Power spectrum baseline, averaged across time points and trials:
%--------------------------------------------------------------------------
tf_pow_baseline = mean(baselinepow,3);
tf_pow_baseline = repmat(tf_pow_baseline,[1 1 size(savetimeidx,2)]);

%--------------------------------------------------------------------------
% Decibel baseline normalization, all time points - averaged across trials:
%--------------------------------------------------------------------------
%tf_pow_db = 10*log10(tf_pow_mean ./ tf_pow_baseline);
%tf_pow_db = 10*log10( bsxfun(@rdivide, tf_pow_mean, tf_pow_baseline) );
tf_pow_db = 10*(log10(tf_pow_mean) - log10(tf_pow_baseline));

%--------------------------------------------------------------------------
% Percent change baseline normalization, all time points - averaged across
% trials:
%--------------------------------------------------------------------------
tf_pow_pc = 100*(tf_pow_mean-tf_pow_baseline)./ tf_pow_baseline;

%--------------------------------------------------------------------------
% Uncomment below to plot results:
%--------------------------------------------------------------------------
% figure
% chanidx = 17;
% contourf(datatimes, frex, mean(mean(squeeze(10*log10(full_pow(chanidx,:,:,:))), 4), 3), 40, 'linecolor', 'none')
% set(gca, 'xlim', [savetime(1) savetime(end)])
% xlabel('Time (ms)'), ylabel('Frequency (Hz)')
% colorbar
%
%
% figure
% chanidx = 52;
% contourf(params.savetime, params.frex, squeeze(tf_pow_db(chanidx,:,:)), 45, 'linecolor', 'none')
% set(gca, 'ydir', 'normal', 'xlim', [-200 1200])
% colormap('jet')
% colorbar
% %caxis([-100 100])
% 
% ch1 = 17;
% ch2 = 52;
% conn = squeeze(squeeze(ispc(ch1, ch2, :, :)));
% 
% for fi=1:length(frex)
% thresh(fi) = median(conn(fi,:)) + std(conn(fi, :));
% tmp = conn(fi,:);
% tmp(tmp < thresh(fi)) = 0;
% connplot(fi,:) = tmp;
% 
% end
% 
% figure
% %chanidx = 60;
% contourf(params.savetime, params.frex, squeeze(conn), 45, 'linecolor', 'none')
% set(gca, 'ydir', 'normal', 'xlim', [-200 1200])
% colormap('jet')
% colorbar
% caxis([0 1])

%%
%--------------------------------------------------------------------------
% Prepare to get cross spectral density for all frequencies, trials and
% time points:
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Preallocate memory for the cross spectral density:
%--------------------------------------------------------------------------
crssspctrm = zeros(dims(1), dims(1), size(frex, 2), size(savetimeidx,2), trials);
blcrssspctrm = zeros(dims(1), dims(1), size(frex, 2), trials);
pli = zeros(dims(1), dims(1), size(frex, 2), size(savetimeidx,2));
ispc = zeros(dims(1), dims(1), size(frex, 2), size(savetimeidx,2));
itpc = zeros(dims(1), size(frex, 2), size(savetimeidx,2));

blpli = zeros(dims(1), dims(1), size(frex, 2));
blispc = zeros(dims(1), dims(1), size(frex, 2));
blitpc = zeros(dims(1), size(frex, 2));

%--------------------------------------------------------------------------
% Because of the complexity of this matrix, the fastest and simplest way to
% calculate the cross spectrum is to pull the channel x channel autospectra
% and use simple array multiplication:
%--------------------------------------------------------------------------
for chani = 1:nbchan
   
    for chanj = chani:nbchan
       
%--------------------------------------------------------------------------
% Update waitbars:
%--------------------------------------------------------------------------
        multiWaitbar( 'Computing Spectral Properties - Channel j', chanj/nbchan );
        multiWaitbar( 'Computing Spectral Properties - Channel i', chani/nbchan );

%--------------------------------------------------------------------------
% Do multiplication of the analytic signal to get the cross spectral
% density matrix for every frequency, time point and trial.
%
% We will get these for both the signal of interest and also the baseline.
%--------------------------------------------------------------------------
        xsd = squeeze(squeeze(bins_analytic_signal(chani,:, :, :) .* conj(bins_analytic_signal(chanj, :, :, :))));
        crssspctrm(chani,chanj, :, :, :) = xsd;

        blxsd = squeeze(squeeze(baseline_analytic_signal(chani,:, :) .* conj(baseline_analytic_signal(chanj, :, :))));
        blcrssspctrm(chani,chanj, :, :) = blxsd;        
       
%--------------------------------------------------------------------------
% Take advantage of the cross spectra - grab a couple of connectivity
% metrics - the phase lag index (PLI) and intersite phase coupling (ISPC):
%--------------------------------------------------------------------------      
        pli(chani,chanj, :, :) = abs(mean(sign(imag(xsd)), 3));
        ispc(chani,chanj, :, :) = abs(mean(exp(1i*angle(xsd)), 3));
       
        blpli(chani,chanj, :) = abs(mean(sign(imag(blxsd)), 2));
        blispc(chani,chanj, :) = abs(mean(exp(1i*angle(blxsd)), 2));

%--------------------------------------------------------------------------
% Because of how this was calculated in a loop instead of multiplication
% of a transpose to get the dot product, let's also populate the lower
% triangles:
%--------------------------------------------------------------------------      
        crssspctrm(chanj,chani, :, :, :) = xsd;
        blcrssspctrm(chanj,chani, :, :) = blxsd;
       
        pli(chanj,chani, :, :) = abs(mean(sign(imag(xsd)), 3));
        ispc(chanj,chani, :, :) = abs(mean(exp(1i*angle(xsd)), 3));
       
        blpli(chanj,chani, :) = abs(mean(sign(imag(blxsd)), 2));
        blispc(chanj,chani, :) = abs(mean(exp(1i*angle(blxsd)), 2));
       
%--------------------------------------------------------------------------
% End for channel j:
%--------------------------------------------------------------------------
    end
   
%--------------------------------------------------------------------------
% Again - take advantage of the cross spectrum - extract inter-trial phase
% clustering (ITPC):
%--------------------------------------------------------------------------
    itpc(chani, :, :) = squeeze(abs(mean(exp(1i*angle(bins_analytic_signal(chani,:, :, :))), 4)));
    blitpc(chani, :) = squeeze(abs(mean(exp(1i*angle(baseline_analytic_signal(chani,:, :, :))), 3)));
   
%--------------------------------------------------------------------------
% End for channel i:
%--------------------------------------------------------------------------
end

%--------------------------------------------------------------------------
% Kill the waitbars:
%--------------------------------------------------------------------------
multiWaitbar('CloseAll');

%--------------------------------------------------------------------------
% Gather output information related to the parameters:
%--------------------------------------------------------------------------
params.datatimes = datatimes;
params.frex = frex;
params.savetime = savetime;
params.basetime = basetime;
params.range_cycles = range_cycles;
params.max_freq = max_freq;
params.min_freq = min_freq;
params.num_frex = size(frex,2);

%--------------------------------------------------------------------------
% Gather all of the spectral power variables:
%--------------------------------------------------------------------------
crss.crssspctrm = crssspctrm;
crss.blcrssspctrm = blcrssspctrm;

%--------------------------------------------------------------------------
% Gather all of the spectral power variables:
%--------------------------------------------------------------------------
%powerspec.tf_pow_trials = tf_pow_trials;
powerspec.tf_pow_mean = tf_pow_mean;
%powerspec.tf_pow_baseline = tf_pow_baseline;
powerspec.tf_pow_db = tf_pow_db;
powerspec.tf_pow_pc = tf_pow_pc;

%--------------------------------------------------------------------------
% Gather all of the spectral power variables:
%--------------------------------------------------------------------------
connmat.pli = pli;
connmat.ispc = ispc;
connmat.blpli = blpli;
connmat.blispc = blispc;

%%