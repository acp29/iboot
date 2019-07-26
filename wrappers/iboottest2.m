%  Function File: iboottest2
%
%  Two-sample bootstrap test for a location parameter
%
%   p = iboottest2(nboot,{bootfun,x,y})
%   p = iboottest2(nboot,{bootfun,x,y},Name,Value)
%   [p,ci] = iboottest2(...)
%
%  Two sample bootstrap test for univariate data. The null hypothesis
%  is that the difference between the bootfun statistic calculated
%  for x and y is equal to zero. The test is two-tailed. Note that
%  this test is only suitable when bootfun is a location parameter.
%  If bootfun is not a location parameter, use boottest2 instead.
%
%  See ibootci documentation for input argument definitions and
%  for Name-Value pairs
%
%  The syntax in this function code is known to be compatible with
%  recent versions of Octave (v3.2.4 on Debian 6 Linux 2.6.32) and
%  Matlab (v7.4.0 on Windows XP).
%
%  iboottest2 v1.0.0.0 (25/07/2019)
%  Author: Andrew Charles Penn
%  https://www.researchgate.net/profile/Andrew_Penn/

function [p,ci] = iboottest2(varargin)

  % Check and process iboottest2 input arguments
  argin = varargin;
  x = argin{2}{2};
  if all(size(x)>1)
    error('data must be a vector')
  end
  if size(x,2)>1
    x = x.';
  end
  y = argin{2}{3};
  if all(size(y)>1)
    error('data must be a vector')
  end
  if size(y,2)>1
    y = y.';
  end
  if numel(argin{2})>3
    error('too many data variables provided')
  end
  nboot = argin{1};
  bootfun = argin{2}{1};
  argin(1:2) = [];

  % Center x and y and calculate vector of pooled residuals
  fx = feval(bootfun,x);
  fy = feval(bootfun,y);
  residuals = cat(1,x-fx,y-fy);

  % Calculate the sample statistic using bootfun
  T0 = fx-fy;

  % Scale and offset the residuals
  residuals = T0 + residuals * 2;

  % Calculate confidence interval using ibootci
  [ci,bootstat,S,calcurve] = ibootci(nboot,{bootfun,residuals},argin{:});

  % Calculate p-value using ibootp
  p = ibootp(0,bootstat,S,calcurve);

end
