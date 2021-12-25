%  Function File: bootperm
%
%  Bootstrap permutation tests
%
%   p = bootperm(y,m)
%   p = bootperm(y,g)
%   p = bootperm(y,...,nboot)
%   p = bootperm(y,...,nboot,bootfun)
%
%   This function provides a bootstrap version of permutation tests [1].
%
%   p = bootperm(y,m) is a 1-sample bootstrap permutation test IF m is
%   a scalar value. m corresponds to the null hypothesis. The default
%   is 0. 
%
%   p = bootperm(y,g) is a k-sample bootstrap permutation test if g is
%   a vector the same size as y with unique numbers used as group labels.
%   If the number of groups (k) is 2, this is a 2-sample bootstrap
%   permutation test. If k > 2, this is a bootstrap permutation test for
%   k groups.
%   
%   p = bootperm(y,...,nboot) sets the number of bootstrap resamples. 
%   The default is 5000.
%
%   p = bootperm(y,...,nboot,bootfun) also sets the statistic calculated
%   from the bootstrap samples. This can be a function handle or string
%   corresponding to the function name. The default is @mean or 'mean'.
%
%   Bibliography:
%   [1] Efron and Tibshirani. Chapter 16 Hypothesis testing with the
%        bootstrap in An introduction to the bootstrap (CRC Press, 1994)
%
%  bootperm v2.0.0.0 (23/12/2021)
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


function p = bootperm(y,vararg,nboot,bootfun)

  % Check and process bootperm input arguments
  if all(size(y)>1) || ~numel(y)>1
    error('y must be a vector')
  end
  if size(y,2)>1
    y = y.'; 
  end
  if all(size(vararg)>1)
    error('the second input argument must be a scalar value (m) or vector (g)')
  end
  if (nargin < 3)
    nboot = 5000;
  end
  if any(size(nboot)>1)
    error('nboot must be scalar. bootperm is not compatible with bootstrap iteration')
  end
  if (nargin < 4)
    bootfun = 'mean';
  end

  % Define function to calculate maximum difference among groups
  if (numel(vararg) > 1)
    % Test for two or more groups 
    g = vararg;
    if size(g,2)>1
      g = g.'; 
    end
    if (numel(g)>1) && (numel(y) ~= numel(g))
      error('y and g vectors must be the same size')
    end
    func = @(y) gfunc(y,g,bootfun);
  else
    % One-sample test
    m = vararg;
    y = cat(1,y,-1*y);
    func = @(y) mfunc(y,bootfun);
  end

  % Perform resampling and calculate bootstrap statistics
  bootstat = bootstrp(nboot,func,y);

  % Calculate p-value
  stat = func(y);
  [cdf,bootstat] = empcdf(bootstat,0);
  p = 1-interp1(bootstat,cdf,stat,'linear');

end

%--------------------------------------------------------------------------

function t = gfunc(Y,g,bootfun)

  % Permutation test statistic for 2 or more groups

  % Get size and of the data vector or matrix
  [m,n] = size(Y);

  % Calculate the number (k) of unique groups
  gk = unique(g);
  k = numel(gk);

  % Calculate maximum difference statistic (t) among the groups
  Z = zeros(k,n);
  for i = 1:k
    Z(i,:) = feval(bootfun,Y(g==gk(i),:));
  end
  Z = sort(Z,1);
  t = Z(k,:)-Z(1,:);

end

%--------------------------------------------------------------------------

function t = mfunc(Y,bootfun)

  % Permutation test statistic for the 1 sample 

  n = size(Y,1)/2;
  Y(n+1:end,:) = [];
  t = abs(feval(bootfun,Y));

end
