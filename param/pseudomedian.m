%  Function file:  [p] = pseudomedian (x, dim)
%
%  If x is a vector, find the univariate pseudo-median of x.
%
%  If x is a matrix, compute the univariate pseudo-median value for
%  each column and return them in a row vector.  If the optional
%  argument dim is given, operate along this dimension. Arrays of
%  more than two dimensions are not currently supported.
%
%  The pseudo-median is the Hodges-Lehmann estimator of centrality,
%  which is robust having a breakdown point of 0.29.
%
%  The syntax in this function code is known to be compatible with
%  recent versions of Octave (v3.2.4 on Debian 6 Linux 2.6.32) and
%  Matlab (v6.5.0 and v7.4.0 on Windows XP).
%
%  pseudomedian v1.0 (last updated: 14/09/2015)
%  Author: Andrew Charles Penn
%  https://www.researchgate.net/profile/Andrew_Penn/


function [p] = pseudomedian (x, dim)

  if nargin<1 || nargin>2
    error('Invalid number of input arguments')
  end

  if ndims(x)>2
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
  if dim == 2
    x = x.';
  end

  % Obtain row dimensions
  s = size(x);
  m = s(1);
  n = s(2);

  % Create xi and xj values with the i <= j restriction enforced
  q = logical(triu(ones(m,m),0));
  i = uint32((1:m)'*ones(1,m));
  xi = x(i(q),:);
  j = uint32(ones(m,1)*(1:m));
  xj = x(j(q),:);

  % Calculate pairwise means (Walsh averages)
  W = (xi+xj)./2;

  % Calculate ordinary median of Walsh averages
  p = median(W);

  % If applicable, switch dimension
  if dim == 2
    p = p.';
  end
