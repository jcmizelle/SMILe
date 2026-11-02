function [cs, ssout] = SMILe_CS_Wavelet(data,segleng,segshift,epleng,band,fs)
%--------------------------------------------------------------------------
% calculates cross-spectrum and coherence based on a wavelet using
% a Hanning window for given frequeny
%
% Input: 
% data:   NxM matrix for N time points and M channels
% segleng:  length of segment (in samples)
% segshift:  shift of segments (in samples)
% epleng: length of epoch (or trial) (in samples)
% band:   frequency band of interest (in Hz)
% fs: sampling frequeny (in Hz)
%
% Output: 
% cs: nchan x nchan cross-spectrum 
% ss: the complex wavelet

%--------------------------------------------------------------------------
% Start calculations:
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Get dimensions for windowing:
%--------------------------------------------------------------------------
npts = (1:segleng)'-segleng/2;

%--------------------------------------------------------------------------
% Create window:
%--------------------------------------------------------------------------
mywin = hanning(segleng);

%--------------------------------------------------------------------------
% Build wavelet:
%--------------------------------------------------------------------------
s1 = cos(npts*band*2*pi/fs).*mywin;
s2 = sin(npts*band*2*pi/fs).*mywin;
ss = s1-sqrt(-1)*s2;
ss = mean(ss,2);
ssout = ss;

[n, nchan] = size(data);
ss = repmat(ss,1,nchan);

%--------------------------------------------------------------------------
% Get parameter definitions:
%--------------------------------------------------------------------------
nep = floor(n/epleng);
nseg = (epleng-segleng)/segshift+1;

%--------------------------------------------------------------------------
% Prepare output:
%--------------------------------------------------------------------------
cs = zeros(nchan,nchan);

kk = 0;
for i = 1:nep
   datalocep = data((i-1)*epleng+1:i*epleng,:);
    for j = 1:nseg
        kk = kk+1;
        datalocseg = datalocep((j-1)*segshift+1:(j-1)*segshift+segleng,:);
        datafreq = transpose(sum(datalocseg.*ss));
        
         if kk == 1
            cs = datafreq*datafreq';
            
         else
            cs = cs+datafreq*datafreq';
            
         end
   
        
    end
    
end

cs = cs/kk;

return;