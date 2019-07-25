%  Function File: boottest2
%
%  Two-sample bootstrap test
%
%   p = iboottest(nboot,bootfun,x,y)
%
%  Two sample bootstrap test for unpaired univariate data. The null
%  hypothesis is that x and y are sampled from the same population.
%  The test is two-tailed.
%
%  The syntax in this function code is known to be compatible with
%  recent versions of Octave (v3.2.4 on Debian 6 Linux 2.6.32) and
%  Matlab (v7.4.0 on Windows XP).
%
%  boottest2 v1.0.0.0 (25/07/2019)
%  Author: Andrew Charles Penn
%  https://www.researchgate.net/profile/Andrew_Penn/

function p = boottest2(nboot,bootfun,x,y)

  % Check and process iboottest2 input arguments
  if any(size(nboot)>1)
    error('boottest2 is not compatible with bootstrap iteration')
  end

  % Calculate confidence interval using ibootci
  z = cat(1,x,y);
  n = numel(x);
  func = @(z) null(bootfun,z,n);
  [ci,bootstat,S,calcurve] = ibootci(nboot,{func,z},'type','per');

  % Calculate p-value using ibootp
  stat = func(z);
  p = ibootp(stat,bootstat,S,calcurve);

end

function t = null(bootfun,z,n)

  % Create difference statistic for the null hypothesis
  x = z(1:n,:);
  y = z(n+1:end,:);
  tx = feval(bootfun,x);
  ty = feval(bootfun,y);
  t = tx-ty;

end
