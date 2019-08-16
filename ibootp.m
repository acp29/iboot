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
%  This functon requires the output bootstrap replicate statistics
%  (bootstat) and the output structure (S) from ibootci. Provision
%  of the calibration curve (calcurve) is optional but highly
%  recommended for accurate P-values.
%
%  P-values obtained from ibootp are consistent with the confidence
%  intervals calculated with ibootci.
%
%  When the p-value is too small to calibrate by double bootstrap,
%  this function automatically resorts to calculating the p-value
%  using the Studentized bootstrap (bootstrap-t) with an additive
%  correction factor to stabilize the variance.
%
%  The syntax in this function code is known to be compatible with
%  recent versions of Octave (v3.2.4 on Debian 6 Linux 2.6.32) and
%  Matlab (v7.4.0 on Windows XP).
%
%  ibootp v1.5.0.0 (16/08/2019)
%  Author: Andrew Charles Penn
%  https://www.researchgate.net/profile/Andrew_Penn/


function p = ibootp(m,bootstat,S,calcurve)

  % Check input arguments
  if nargin < 2
    error('ibootp requires atleast 2 input arguments')
  end
  if nargin < 3
    S.z0 = 0;
    S.a = 0;
  end

  % Calculate number of bootstrap replicates
  if iscell(bootstat)
    T1 = bootstat{1};
  else
    T1 = bootstat;
  end
  B = numel(T1);

  % Find P-value from the cdf of the bootstrap distribution by linear interpolation
  [cdf,t1] = empcdf(T1,0);
  p = 1-interp1(t1,cdf,m,'linear','extrap');

  % BCa correction to P-value if applicable
  z1 = norminv(p);
  p = normcdf(S.z0+((S.z0+z1)/(1-S.a*(S.z0+z1))));

  % Convert to equal-sided two-tailed test
  p = 2*min(p,1-p);

  % Check if first bootstrap replicate sample set is large enough
  if (1/min(p,1-p)) > (0.5*B) || isnan(p)
    warning(sprintf(['P value too small for this bootstrap distribution. \n'...
            'Try increasing the number of first bootstrap replicate samples in ibootci.']));
    if isnan(p)
      p = 0;
    end
  end

  % Calibration of P-value if applicable
  if nargin > 3
    C = S.nboot(2);
    if C > 0
      if (1/min(p,1-p)) < (0.5*C)

        % Use same calibration of p value as used for confidence intervals
        calcurve(1,:)=[];calcurve(end,:)=[];
        p = 1 - interp1(calcurve(:,1),calcurve(:,2),1-p,'linear');

      else

        % Use bootstrap-t method with variance stabilization for small samples
        % Polansky (2000) Can J Stat. 28(3):501-516
        se = std(bootstat{1});
        a = S.n^(-3/2)*se;  % additive correction factor

        % Calculate Studentized statistics
        T = (bootstat{1} - S.stat)./(a + std(bootstat{2}));
        t = (S.stat - m)/(a + se);

        % Calculate p value from empirical distribution of the
        % Studentized bootstrap statistics
        [cdf,T] = empcdf(T,0);
        p = 1-interp1(T,cdf,t,'linear','extrap');
        p = 2*min(p,1-p);

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
