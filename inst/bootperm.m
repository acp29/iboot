%  Function File: bootperm
%
%  Bootstrap permutation tests
%
%  p = bootperm(y,m)
%  p = bootperm(y,g)
%  p = bootperm(y,...,nboot)
%  p = bootperm(y,...,nboot,bootfun)
%  p = bootperm(y,...,nboot,bootfun,paropt)
%
%  This function provides a bootstrap version of permutation tests for
%  univariate data [1].
%
%  p = bootperm(y,m) is a 1-sample bootstrap permutation test IF m is
%  a scalar value. m corresponds to the null hypothesis. The default
%  is 0. 
%
%  p = bootperm(y,g) is a k-sample bootstrap permutation test if g is
%  a vector the same size as y with unique numbers used as group labels.
%  If the number of groups (k) is 2, this is a 2-sample bootstrap
%  permutation test. If k > 2, this is a bootstrap permutation test for
%  k groups.
%   
%  p = bootperm(y,...,nboot) sets the number of bootstrap resamples. 
%  The default is 5000.
%
%  p = bootperm(y,...,nboot,bootfun) also sets the statistic calculated
%  from the bootstrap samples. This can be a function handle or string
%  corresponding to the function name. The default is @mean or 'mean'.
%
%  p = bootperm(y,...,nboot,bootfun,paropt) specifies options that govern 
%  if and how to perform bootstrap iterations using multiple processors 
%  (if the Parallel Computing Toolbox or Octave Forge package is available).
%  This argument is a structure with the following recognised fields:
%
%   'UseParallel' - If true, compute bootstrap iterations in parallel.
%                   Default is false for serial computation. In MATLAB,
%                   the default is true if a parallel pool has already
%                   been started.
%   'nproc'       - The number of processors to use by Octave. Default
%                   is the number of available processors. If you choose
%                   In Matlab, nproc is ignored and the number of parallel
%                   workers should be predefined beforehand by starting
%                   a parallel pool, else it will use the preferred number
%                   of workers.
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


function p = bootperm(y,vararg,nboot,bootfun,paropt)

  % Check and process bootperm input arguments
  if all(size(y)>1) || ~numel(y)>1
    error('y must be a vector')
  end
  if size(y,2)>1
    y = y.'; 
  end
  if (nargin < 2)
    vararg = 0;
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
  if (nargin < 5)
    paropt = struct;
    paropt.UseParallel = false;
    % Initialise nproc if it doesn't exist
    if ~exist('nproc','builtin')
      nproc = 1;
    end
    paropt.nproc = nproc;
  end

  % Define function to calculate maximum difference among groups
  if (numel(vararg) > 1)
    % Test for two or more groups
    % Data is exchangeable across all the groups labelled in g
    g = vararg;
    if size(g,2)>1
      g = g.'; 
    end
    if (numel(g)>1) && (numel(y) ~= numel(g))
      error('y and g must be vectors the same size')
    end
    func = @(y) gfunc(y,g,bootfun);
  else
    % One-sample test 
    % Data is exchangeable with symmetry above and below m
    m = vararg;
    y = cat(1,y-m,m-y);
    func = @(y) mfunc(y,bootfun);
  end

  % Perform resampling and calculate bootstrap statistics
  [~,bootstat] = ibootci(nboot,{func,y},'Options',paropt);

  % Calculate p-value
  stat = func(y);
  p = sum(bootstat>stat)/nboot;
  if (p == 0)
    warning('p-value too small to calculate. Try increasing nboot.')
  end

end

%--------------------------------------------------------------------------

function t = gfunc(Y,g,bootfun)

  % Permutation test statistic for 2 or more groups

  % Get size and of the data vector or matrix
  [m,n] = size(Y);

  % Calculate the number (k) of unique groups
  gk = unique(g);
  k = numel(gk);

  % Calculate bootfun on the data from each group
  % (which is simpler and more intuitive than calculating F or MSE)
  Z = zeros(k,n);
  for i = 1:k
    Z(i,:) = feval(bootfun,Y(g==gk(i),:));
  end

  % Calculate maximum difference statistic (t) among the groups
  Z = sort(Z,1);
  t = Z(k,:)-Z(1,:); % sign always positive

end

%--------------------------------------------------------------------------

function t = mfunc(Y,bootfun)

  % Permutation test statistic for the 1 sample case
  
  % Make sample the same size as the original data 
  n = size(Y,1)/2;
  Y(n+1:end,:) = [];

  % Calculate absolute value of result of bootfun on the resample(s)
  t = abs(feval(bootfun,Y)); % sign always positive

end
