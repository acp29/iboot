%  Function File: iboottest2
%
%  Two-sample bootstrap test
%
%   p = iboottest2(nboot,{bootfun,x,y})
%   p = iboottest2(nboot,{bootfun,x,y},Name,Value)
%   [p,ci] = iboottest2(...)
%   [p,ci,S] = iboottest2(...)
%
%  Two sample bootstrap test for univariate data. The null hypothesis
%  is that the difference between the bootfun statistic calculated
%  for x and y is equal to zero. The test is two-tailed.
%
%  See ibootci documentation for input argument definitions and
%  for Name-Value pairs
%
%  The syntax in this function code is known to be compatible with
%  recent versions of Octave (v3.2.4 on Debian 6 Linux 2.6.32) and
%  Matlab (v7.4.0 on Windows XP).
%
%  iboottest2 v1.1.0.0 (25/07/2019)
%  Author: Andrew Charles Penn
%  https://www.researchgate.net/profile/Andrew_Penn/

function [p,ci,S] = iboottest2(argin1,argin2,varargin)

  % Check and process iboottest2 input arguments
  nboot = argin1;
  bootfun = argin2{1};
  argin = varargin;
  x = argin2{2};
  if all(size(x)>1)
    error('data must be a vector')
  end
  if size(x,2)>1
    x = x.';
  end
  y = argin2{3};
  if all(size(y)>1)
    error('data must be a vector')
  end
  if size(y,2)>1
    y = y.';
  end
  if numel(argin2)>3
    error('too many data variables provided')
  end

  % Perform independent resampling from x and y
  state = warning;
  warning off;
  [~,bootstatX,SX] = ibootci(nboot,{bootfun,x},argin{:});
  [~,bootstatY,SY] = ibootci(nboot,{bootfun,y},argin{:});
  warning(state);

  % Calculate differences between bootfun evaluated for bootstrap
  % sample sets derived from x and y
  if iscell(bootstatX)
    bootstatZ = cell(2,1);
    bootstatZ{1} = bootstatX{1} - bootstatY{1};
    bootstatZ{2} = bootstatX{2} - bootstatY{2};
  else
    bootstatZ = bootstatX - bootstatY;
  end

  % Create template settings structure and calculate the sample test statistic
  S = SX;
  T0 = SX.stat - SY.stat;
  S.stat = T0;

  % Calculate confidence interval using ibootci
  [ci,bootstat,S,calcurve] = ibootci(bootstatZ, S);

  % Calculate p-value using ibootp
  p = ibootp(0,bootstat,S,calcurve);

end


