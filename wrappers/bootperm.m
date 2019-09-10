%  Function File: bootperm
%
%  One-sample bootstrap (permutation) test
%
%   p = bootperm(nboot,bootfun,x,m)
%   p = bootperm(nboot,bootfun,x,m,clusters)
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
%  The optional input argument, clusters must be provided as a cell
%  array. The first cell should contain cluster definitions for x 
%  and the second cell should contain cluster definitions for y. 
%  An empty cell signifies that no clusters will be used in the 
%  bootstrap for that sample. See ibootci documentation for how to 
%  construct cluster definitions.
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


function p = bootperm1(nboot,bootfun,x,m,clusters)

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
  global J
  if nargin > 4
    if ~isempty(clusters)
      if ~all(size(clusters) == size(x))
        error('dimensions of the data and clusters must be the same')
      end
      % Concatenate clusters
      clusters = repmat(clusters,2,1);
      % Sort clusters
      [clusters,I] = sort(clusters);
      [~,J] = sort(I);
    end
  else
    clusters = [];
    I = (1:2*n)';
    J = I;
  end 

  % Prepare joint distribution (and sort to match clusters if necessary) 
  z = cat(1,x,-1*x);
  z = z(I);
  
  % Define function to test the null hypothesis
  func = @(z) null(bootfun,z,n);

  % Use ibootci to create bootstrap statistics
  [~,bootstat] = ibootci(nboot,{func,z},'type','per','Clusters',clusters);

  % Calculate p-value using ibootp
  stat = func(z);
  p = ibootp(stat,bootstat);
  
end
        
%--------------------------------------------------------------------------

function t = null(bootfun,z,n)

  % Resort resampled data
  global J
  z = z(J);
  
  % Calculate statistic for the null hypothesis
  x = z(1:n,:);
  t = feval(bootfun,x);

end
