%  Function file: smoothmedian
%
%  [p] = smoothmedian (x)
%  [p] = smoothmedian (x, dim)
%  [p] = smoothmedian (x, dim, Tol)
%
%  If x is a vector, find the univariate smoothed median of x.
%
%  If x is a matrix, compute the univariate smoothed median value
%  for each column and return them in a row vector.  If the optional
%  argument dim is given, operate along this dimension. Arrays of
%  more than two dimensions are not currently supported.
%
%  Smoothing the median is achieved by minimizing the following
%  objective function:
%
%        S (p) = sum {(x(i) - p).^2 + (x(j) - p).^2} .^ 0.5
%               i < j
%
%  With the ordinary median as the initial value of p, this function
%  minimizes the objective function by finding the root of the first
%  derivative using a Newton-Bisection hybrid algorithm. By default,
%  the tolerance (Tol) for the first derivative is set at single
%  machine precision.
%
%  The smoothing works by slightly reducing the breakdown point
%  of the median. Bootstrap confidence intervals using the smoothed
%  median have good coverage for the ordinary median of the
%  population distribution and can be used to obtain second order
%  accurate intervals with Studentized bootstrap and calibrated
%  percentile bootstrap methods [1]. These bootstrap methods are
%  available in bootci (from the Statistics and Machine Learning
%  Toolbox) and ibootci (from Matlab Central File Exchange)
%  respectively. When the population distribution is thought to be
%  strongly skewed, coverage errors can be reduced by improving
%  symmetry through appropriate data transformation. Unlike kernel-
%  based smoothing approaches, bootstrapping smoothmedian does not
%  require explicit choice of a smoothing parameter or a probability
%  density function. The algorithm used is suitable for small-to-
%  medium sample sizes.
%
%  Bibliography:
%  [1] Brown, Hall and Young (2001) The smoothed median and the
%       bootstrap. Biometrika 88(2):519-534
%
%  The syntax in this function code is known to be compatible with
%  recent versions of Octave (v3.2.4 on Debian 6 Linux 2.6.32) and
%  Matlab (v6.5.0 and v7.4.0 on Windows XP).
%
%  smoothmedian v1.4.4 (20/08/2019)
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


function p = smoothmedian(x,dim,Tol)

  if nargin<1 || nargin>3
    error('Invalid number of input arguments')
  end

  if ~ismatrix(x)
    error('arrays of more than 2 dimensions are not supported')
  end

  if nargin<2
    if size(x,2)==1
      dim = 1;
    elseif size(x,1)==1
      dim = 2;
    else
      dim = 1;
    end
  end

  if dim<1 || dim>2
    error('dim must be a valid dimension');
  end

  % If applicable, switch dimension
  if dim > 1
    x = x.';
  end

  if nargin<3
    Tol = 1.1921e-07;
  end

  if ~isa(x,'double')
    error('the x variable must be double precision')
  end

  % Obtain data dimensions
  s = size(x);
  m = s(1);
  n = s(2);
  l = m*(m-1)/2;

  % Find column indices where smoothing is not possible
  no = any(isnan(x))|any(isinf(x));

  % Variable transformation to normalize stopping criteria
  % Centre the data on the median and divide by the midrange
  % Regularization is achieved by adding a constant to the denominator
  centre = nanmedian(x,1);
  midrange = (max(x,[],1)-min(x,[],1))/2;
  x = (x-centre(ones(1,m),:))./(1+midrange(ones(1,m),:));

  % Obtain m(m-1)/2 pairs from the Cartesian product of each column of
  % x with itself by enforcing the restriction i < j on xi and xj
  q = logical(triu(ones(m,m),1));
  i = uint32((1:m)'*ones(1,m));
  xi = x(i(q),:);
  i = []; %#ok<NASGU> Reduce memory usage. Faster than using clear.
  j = uint32(ones(m,1)*(1:m));
  xj = x(j(q),:);
  j = []; %#ok<NASGU> Reduce memory usage. Faster than using clear.
  q = []; %#ok<NASGU> Reduce memory usage. Faster than using clear.

  % Offset xi and xj by +/- h, where h is small enough to avoid
  % encountering stationary points when calculating function
  % derivatives for the Newton steps.
  h = sqrt(0.5*eps^(2/3));
  xi = xi+h;
  xj = xj-h;
  z = (xi-xj).^2;
  y = xi+xj;

  %% Nonlinear root finding by Newton-Bisection hybrid algorithm
  % Set initial bracket bounds
  a = min(x,[],1);
  b = max(x,[],1);
  % Set starting value as the median
  c = zeros(1,n);
  p = c;
  % Initialize and define settings
  idx = 1:n;
  MaxIter = 5000;
  % Start iterations
  for Iter = 1:MaxIter
    % Compute derivatives
    temp = ones(l,1)*c;
    D = (xi-temp).^2+(xj-temp).^2;
    R = sqrt(D);
    T = sum((2*temp-y)./R,1);
    if Iter==1
      T(no) = 0;
    end
    temp = []; %#ok<NASGU> Reduce memory usage. Faster than using clear.
    U = sum(z.*R./D.^2,1);
    D = []; %#ok<NASGU> Reduce memory usage. Faster than using clear.
    R = []; %#ok<NASGU> Reduce memory usage. Faster than using clear.
    % Update bracket bounds
    a(T<-Tol) = c(T<-Tol);
    b(T>+Tol) = c(T>+Tol);
    % Evaluate convergence
    cvg = abs(T)<Tol;
    if any(cvg)
      % Export converged parameters
      p(idx(cvg)) = c(cvg);
      % Avoid excess computations in following iterations
      idx(cvg) = [];
      xi(:,cvg) = [];
      xj(:,cvg) = [];
      z(:,cvg) = [];
      y(:,cvg) = [];
      a(cvg) = [];
      b(cvg) = [];
      c(cvg) = [];
      T(cvg) = [];
      U(cvg) = [];
    elseif all(cvg)
      % Export converged parameters and break from iterations
      p(idx(cvg)) = c(cvg);
      break
    end
    % Compute Newton step (fast quadratic convergence but unreliable)
    W = T./U;
    d = c-W;
    T = []; %#ok<NASGU> Reduce memory usage. Faster than using clear.
    U = []; %#ok<NASGU> Reduce memory usage. Faster than using clear.
    W = []; %#ok<NASGU> Reduce memory usage. Faster than using clear.
    % Compute Bisection step (slow linear convergence but very safe)
    c = 0.5*(a+b);
    % Prefer Newton step if it is within brackets
    nwt = d>a & d<b;
    c(nwt) = d(nwt);
    d = [];  %#ok<NASGU> Reduce memory usage. Faster than using clear.
    nwt =[]; %#ok<NASGU> Reduce memory usage. Faster than using clear.
  end
  if Iter==MaxIter
    fprintf('Warning: root finding failed to reach the specified tolerance.\n');
  end

  % Backtransform the smoothed median value(s)
  p = p.*(1+midrange)+centre;

  % Assign ordinary median where smoothing is not possible
  p(no) = centre(no);

  % If applicable, switch dimension
  if dim > 1
    p = p.';
  else
    p = p;
  end
