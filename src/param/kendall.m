%  Function file: kendall
%
%  [tau] = kendall (x,y)
%  [tau] = kendall (x,y,dim)
%  [tau,p] = kendall (...)
%
%  This function returns Kendall's tau rank correlation coefficient for x
%  and y [1]. If x and y are column vectors, a scalar value is returned.
%  If x and y are matrices, the statistic is computed for matching columns
%  and the result is returned as a row vector. If the optional argument
%  dim is given, operate along this dimension. This function computes the
%  tau-b statistic adjusted for ties [2,3] using a vectorized form of the
%  follow algorithm:
%
%      tau  = S ./ D
%
%    Where,
%
%      S  =  sum {sgn(x(i) - x(j)) .* sgn(y(i) - y(j))}
%           i < j
%
%      D  = (sum 1{|x(i) - x(j)| > 0} .* sum 1{|y(i) - y(j)| > 0}).^0.5
%           i < j                       i < j
%
%      1 is the indicator function
%
%  The resulting distribution free statistic is a value between -1 and
%  +1 that reflects the nature and strength of dependence between two
%  variables. Although any relationship between the two variables is
%  not presumed linear, the statistic does assume that they vary
%  monotonically.
%
%  The statistic can be also be represented as:
%
%    tau   =       ((# of concordant pairs) - (# of discordant pairs))
%            ----------------------------------------------------------
%            sqrt ((# of untied pairs in x) * (# of untied pairs in y))
%
%
%  The algorithm used is suitable for small-to-mediam sample sizes.
%
%  This function can also compute two-tailed p-values for Kendall's 
%  tau. The p-value is calculated for small samples (< 9) by exact 
%  permutation. For larger samples the p-value is approximated from 
%  a sampling distribution estimated using 50,000 random permutations.
%
%  The syntax in this function code is known to be compatible with
%  recent versions of Octave (v3.2.4 on Debian 6 Linux 2.6.32) and
%  Matlab (v6.5.0 and v7.4.0 on Windows XP).
%
%  Bibliography:
%  [1] Kendall, M. (1938) A New Measure of Rank Correlation. Biometrika
%       30 (1-2): 81-89.
%  [2] Kendall, M. G. (1945). The treatment of ties in ranking problems.
%       Biornetrika 33, 239-51
%  [3] Kendall, M. G. (1947). The variance of tau when both rankings
%       contain ties. Biometrika 34, 297-8. 
%
%  kendall v1.0 (last updated: 19/08/2013)
%  Author: Andrew Charles Penn
%  https://www.researchgate.net/profile/Andrew_Penn/
%
%  Copyright 2010 Andrew Charles Penn
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


function [tau,p] = kendall(x,y,dim)

  if nargin<2
    error('Too few input arguments');
  end

  if nargout>2
    error('Too many output arguments');
  end

  if nargin<3
    dim = 1;
  end

  if dim<1 || dim>2
    error('dim must be a valid dimension');
  end

  if ~isa(x,'numeric')
    error('The x variable must be numeric')
  end

  if ~isa(y,'numeric')
    error('The y variable must be numeric')
  end

  if ~all(size(x)==size(y))
    error('The dimensions of x and y are not consistent')
  end

  if max(size(x))==1
    error('x and y must be vectors or matrices')
  end

  % If applicable, switch dimension
  if dim==2
    x = x';
    y = y';
  end
  
  % Create xi, xj, yi and yj with the i < j restriction enforced
  m = size(x,1);
  q = logical(triu(ones(m,m),1));
  i = uint32((1:m)'*ones(1,m));
  xi = x(i(q),:);
  yi = y(i(q),:);
  i = []; %#ok<NASGU> Reduce memory usage. Faster than using clear.
  j = uint32(ones(m,1)*(1:m));
  xj = x(j(q),:);
  yj = y(j(q),:);
  j = []; %#ok<NASGU> Reduce memory usage. Faster than using clear.
  q = []; %#ok<NASGU> Reduce memory usage. Faster than using clear.
  
  % Calculate differences between x and y values at positions i and j 
  xdiff = xi-xj;
  ydiff = yi-yj;
  xi = []; %#ok<NASGU> Reduce memory usage. Faster than using clear.
  xj = []; %#ok<NASGU> Reduce memory usage. Faster than using clear.
  yi = []; %#ok<NASGU> Reduce memory usage. Faster than using clear.
  yj = []; %#ok<NASGU> Reduce memory usage. Faster than using clear.
  
  % Calculate Kendall's tau-b correlation coefficient adjusted for ties 
  S = sum(sign(xdiff).*sign(ydiff),1);
  D = sqrt(sum(abs(xdiff)>0,1).*sum(abs(ydiff)>0,1));
  tau = S./D;
  
  % Compute two-tailed p-value by permutation
  % Although permutation is slow, it does not use the Normal approximation
  if nargout>1
    if m <= 8
      % For small size samples compute p-value using exact permutation
      X = perms(x)';
      N = size(X,2);
      H0 = kendall(X,y*ones(1,N));
    else
      % For larger samples, estimate the sampling distribution under the
      % null hypothesis using 50,000 random permutations and compute an
      % approximate p-value
      N = 5e+05;
      H0 = zeros(N,1);
      for k = 1:N
        H0(k) = kendall(x(randperm(m)),y);
      end
    end
    p = sum(abs(H0)>=abs(tau))/N;
  end
             
  % If applicable, switch back dimension
  if dim==2
    tau = tau';
  end