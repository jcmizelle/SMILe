function [CI, CONNsd] = SMILe_CIcalc(connmat)
%--------------------------------------------------------------------------
%
% Compute the 95% CI for connectivity data. 
%
% This requires the statistics toolbox for tinv.
%
% Input should be a single nchan x nchan matrix of connectivity. Put this
% in a loop (one matrix at a time) to expand.
%
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Rearrange data to lose one of the nchan dimensions - which is redundant
% and unfairly contributes to CI calculations. 
%--------------------------------------------------------------------------
CIdat = connmat(find(triu(connmat)));


%--------------------------------------------------------------------------
% Get data dimensions. 
%--------------------------------------------------------------------------
CIdims = size(CIdat);
CInumel = prod(CIdims);

%--------------------------------------------------------------------------
% Combine all values, skipping the redundant dimension:
%--------------------------------------------------------------------------
CIdat = reshape(CIdat, 1, CInumel)';  

%--------------------------------------------------------------------------
% Compute standard error of the mean:
%--------------------------------------------------------------------------
CIdatsd = std(CIdat);
%--------------------------------------------------------------------------
% Compute standard error of the mean:
%--------------------------------------------------------------------------
SEM = CIdatsd/sqrt(length(CIdat));

%--------------------------------------------------------------------------
% Estimate t-score:
%--------------------------------------------------------------------------
tscore = tinv([0.025  0.975], length(CIdat)-1);

%--------------------------------------------------------------------------
% Calculate confidence interfals:
%--------------------------------------------------------------------------
CI = mean(CIdat) + tscore * SEM;

%--------------------------------------------------------------------------
% Also export mean +/- 1, 2 and 3 sd:
%--------------------------------------------------------------------------
CONNsd = [mean(CIdat) - CIdatsd mean(CIdat) + CIdatsd;
          mean(CIdat) - CIdatsd*2 mean(CIdat) + CIdatsd*2;
          mean(CIdat) - CIdatsd*3 mean(CIdat) + CIdatsd*3];
        
%--------------------------------------------------------------------------
% End for parent function:
%--------------------------------------------------------------------------
end