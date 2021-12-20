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
%  The test is two-tailed and is related to a one-sample permutation
%  test. Note that this function resamples the rows of the data x.
%
%  bootperm v1.1.2.0 (01/10/2021)
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

  % Get size of data
  n = size(x,1);

  % Evaluate optional input arguments for nested data structures
  if nargin > 3
    if isempty(m)
      m = 0;
    end
    if ~isa(m,'numeric') || numel(m)~=1
      error('m must be a numeric scalar value');
    end
    x = x - m;
  end

  % Prepare joint distribution and define function to test the null hypothesis
  z = cat(1,x,-1*x);
  func = @(z) null(bootfun,z,n);

  % Use ibootci to create bootstrap statistics
  [ci,bootstat] = ibootci(nboot,{func,z},'type','per');

  % Calculate p-value using ibootp
  stat = func(z);
  p = ibootp(stat,bootstat);

end

%--------------------------------------------------------------------------

function t = null(bootfun,z,n)

  % Calculate statistic for the null hypothesis
  x = z(1:n);
  t = feval(bootfun,x);

end
