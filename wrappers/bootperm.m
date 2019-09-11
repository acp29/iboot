%  Function File: bootperm
%
%  One-sample bootstrap (permutation) test
%
%   p = bootperm(nboot,bootfun,x,m)
%   p = bootperm(nboot,bootfun,x,m,strata)
%
%  One sample bootstrap test for univariate data. The null hypothesis
%  is that x is sampled from a population with location m calculated
%  by bootfun. This test performs well without bootstrap iteration
%  and so iteration is not implemented but it can only compute a P
%  value (not a confidence interval).
%
%  The test is two-tailed and is related to a one-sample permutation 
%  test.
%
%  See ibootci documentation for information about the optional input 
%  argument strata.
%
%  The syntax in this function code is known to be compatible with
%  recent versions of Octave (v3.2.4 on Debian 6 Linux 2.6.32) and
%  Matlab (v7.4.0 on Windows XP).
%
%  bootperm v1.1.0.0 (10/09/2019)
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


function p = bootperm(nboot,bootfun,x,m,strata)

  % Check and process bootperm input arguments
  if any(size(nboot)>1)
    error('bootperm is not compatible with bootstrap iteration')
  end
  
  % Get size of data
  n = numel(x);  
  
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
  if nargin > 4
    if ~isempty(strata)
      if ~all(size(strata) == size(x))
        error('dimensions of the data and strata must be the same')
      end
      % Sort strata definitions and data
      [strata,I] = sort(strata);
      x = x(I);
    end
  else
    strata = [];
  end 

  % Define function to test the null hypothesis
  func = @(x) null(bootfun,x,n);

  % Use ibootci to create bootstrap statistics
  [~,bootstat] = ibootci(nboot,{func,x},'type','per','Strata',strata);
  
  % Calculate p-value using ibootp
  stat = func(x);
  p = ibootp(stat,bootstat);
  
end
        
%--------------------------------------------------------------------------

function t = null(bootfun,x,n)

  % Calculate statistic for the null hypothesis
  z = bsxfun (@times, x, 2*(randi(2,n,1)-1.5));
  t = feval(bootfun,z);

end
