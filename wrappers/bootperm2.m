%  Function File: bootperm2
%
%  Two-sample bootstrap (permutation) test
%
%   p = bootperm2(nboot,bootfun,x,y)
%
%  Two sample bootstrap test for unpaired univariate data. The null
%  hypothesis is that x and y are sampled from the same population.
%  The test is two-tailed and is essentially a bootstrap version of
%  a permutation test. This test performs well without bootstrap
%  iteration and so iteration is not implemented but it can only
%  compute a P value (not a confidence interval).
%
%  The syntax in this function code is known to be compatible with
%  recent versions of Octave (v3.2.4 on Debian 6 Linux 2.6.32) and
%  Matlab (v7.4.0 on Windows XP).
%
%  bootperm2 v1.0.2.0 (15/08/2019)
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


function p = bootperm2(nboot,bootfun,x,y)

  % Check and process ibootperm2 input arguments
  if any(size(nboot)>1)
    error('bootperm2 is not compatible with bootstrap iteration')
  end

  % Prepare joint distribution and define function to test the null hypothesis
  z = cat(1,x,y);
  n = numel(x);
  func = @(z) null(bootfun,z,n);

  % Use ibootci to create bootstrap statistics
  state = warning;
  warning off;
  [~,bootstat] = ibootci(nboot,{func,z},'type','per');
  warning(state);

  % Calculate p-value using ibootp
  stat = func(z);
  p = ibootp(stat,bootstat);

end

%--------------------------------------------------------------------------

function t = null(bootfun,z,n)

  % Create difference statistic for the null hypothesis
  x = z(1:n,:);
  y = z(n+1:end,:);
  tx = feval(bootfun,x);
  ty = feval(bootfun,y);
  t = tx-ty;

end
