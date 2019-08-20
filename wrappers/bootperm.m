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
%  bootperm v1.0.1.0 (15/08/2019)
%  Author: Andrew Charles Penn
%  https://www.researchgate.net/profile/Andrew_Penn/
%
%  Copyright 2019 Andrew Charles Penn
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <http://www.gnu.org/licenses/>.


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
