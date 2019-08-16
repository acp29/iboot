%  Function File: bootperm
%
%  One-sample bootstrap (permutation) test
%
%   p = bootperm(nboot,bootfun,x,m)
%
%  One sample bootstrap test for univariate data. The null hypothesis
%  is that x is sampled from a population with location m calculated
%  by bootfun. This test performs well without bootstrap iteration
%  and so iteration is not implemented but it can only compute a P
%  value (not a confidence interval).
%
%  The test is two-tailed and is essentially a bootstrap version of
%  a one-sample permutation test.
%
%  The syntax in this function code is known to be compatible with
%  recent versions of Octave (v3.2.4 on Debian 6 Linux 2.6.32) and
%  Matlab (v7.4.0 on Windows XP).
%
%  bootperm v1.0.0.0 (15/08/2019)
%  Author: Andrew Charles Penn
%  https://www.researchgate.net/profile/Andrew_Penn/

function p = bootperm(nboot,bootfun,x,m)

  % Check and process bootperm input arguments
  if any(size(nboot)>1)
    error('bootperm is not compatible with bootstrap iteration')
  end
  if nargin < 4
    m = 0;
  else
    if ~isa(m,'numeric') || numel(m)~=1
      error('m must be a numeric scalar value');
    end
    x = x - m;
  end

  % Prepare vector of signs
  n = numel(x);
  z = ones(n,1);
  k = fix(n/2);
  z(1:k) = -1;

  % If n is odd, create weights to sample signs equally
  isodd = sum(z);
  w = ones(n,1);
  if isodd
    w(1:k) = 0.5/sum(z==-1);
    w(k+1:n) = 0.5/sum(z==1);
  end

  % Use ibootci to create bootstrap indices
  state = warning;
  warning off;
  [ci,bootstat,S,calcurve,bootidx] = ibootci(nboot,{bootfun,z},'Weights',w);
  warning(state);

  % Apply bootstrapped signs to data vector x and calculate bootfun
  Z = z(bootidx);
  bootstat = feval(bootfun,bsxfun(@times,x,Z));

  % Calculate p-value using ibootp
  stat = bootfun(x);
  p = ibootp(stat,bootstat);

end
