%  Function File: ibootp
%
%  Bootstrap P-value: One sample test
%
%   p = ibootp(m,bootstat,S,calcurve)
%
%  Performs a two-tailed bootstrap test for the hypothesis that
%  the data used to generate the bootstrap statistics bootstat
%  comes from a distribution with statistic m. m must be a scalar.
%
%  Requires the output bootstrap replicate statistics (bootstat)
%  and the output structure (S) from ibootci. Provision of the
%  calibration curve (calcurve) is optional but recommended for
%  accurate P-values.
%
%  Unlike for confidence intervals, large numbers of bootstrap
%  replicate samples for both the first and second bootstrap
%  replicate sets are required to generate calibrated bootstrap
%  P-values. The smaller the P-value, the greater the required
%  replicate sample sizes. Warnings are provided when ibootci
%  needs to be run again with larger replicate sample sets.
%
%  P-values obtained from ibootp are consistent with the
%  confidence intervals obtained with ibootci.
%
%  The syntax in this function code is known to be compatible with
%  recent versions of Octave (v3.2.4 on Debian 6 Linux 2.6.32) and
%  Matlab (v7.4.0 on Windows XP).
%
%  ibootp v1.4.2.0 (25/07/2019)
%  Author: Andrew Charles Penn
%  https://www.researchgate.net/profile/Andrew_Penn/


function p = ibootp(m,bootstat,S,calcurve)

  if nargin < 3
    error('ibootp requires atleast 3 input arguments')
  end

  % Calculate number of bootstrap replicates
  B = numel(bootstat);

  % Find P-value from the cdf of the bootstrap distribution by linear interpolation
  [cdf,t1] = empcdf(bootstat,0);
  p = 1-interp1(t1,cdf,m,'linear','extrap');

  % BCa correction to P-value if applicable
  z1 = norminv(p);
  p = normcdf(S.z0+((S.z0+z1)/(1-S.a*(S.z0+z1))));

  % Convert to equal-sided two-tailed test
  p = 2*min(p,1-p);

  % Check if first bootstrap replicate sample set is large enough
  if (1/min(p,1-p)) > (0.5*B) || isnan(p)
    warning(['P value too small for this bootstrap distribution. \n'...
            'Try increasing the number of first bootstrap replicate samples in ibootci.']);
  end

  % Calibration of P-value if applicable
  if nargin > 3
    C = S.nboot(2);
    if C > 0
      calcurve(1,:)=[];calcurve(end,:)=[];
      if (1/min(p,1-p)) < (0.5*C)
        p = 1 - interp1(calcurve(:,1),calcurve(:,2),1-p,'linear');
      else
        warning(['P value is too small for calibration so the result will be unreliable.\n'...
                'Try increasing the number of second bootstrap replicate samples in ibootci.']);
      end
    end
  end

end

%--------------------------------------------------------------------------

function [F,x] = empcdf (y,c)

  % Calculate empirical cumulative distribution function of y
  %
  % Set c to:
  %  1 to have a complete distribution with F ranging from 0 to 1
  %  0 to avoid duplicate values in x
  %
  % Unlike ecdf, empcdf uses a denominator of N+1

  % Check input argument
  if ~isa(y,'numeric')
    error('y must be numeric')
  end
  if all(size(y)>1)
    error('y must be a vector')
  end
  if size(y,2)>1
    y = y.';
  end

  % Create empirical CDF
  x = sort(y);
  N = sum(~isnan(y));
  [x,F] = unique(x,'rows','last');
  F = F/(N+1);

  % Apply option to complete the CDF
  if c > 0
    x = [x(1);x;x(end)];
    F = [0;F;1];
  end

end
