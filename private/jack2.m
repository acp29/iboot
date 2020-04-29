function [SE, T, U] = jack (x, func, paropt)

  % Private function file required for bootci

  % Ordinary Jackknife

  if nargin < 2
    error('Invalid number of input arguments');
  elseif nargin < 3
    paropt = struct;
    paropt.UseParallel = false;
    paropt.nproc = 1;
  end

  if nargout > 3
    error('Invalid number of output arguments');
  end

  % Perform 'leave one out' procedure and calculate the variance(s)
  % of the test statistic.
  nvar = size(x, 2);
  m = size(x{1}, 1);
  ridx = diag(ones(m, 1));
  j = (1:m)';
  M = cell(1, nvar);
  for v = 1:nvar
    M{v} = x{v}(j(:, ones(m, 1)), :);
    M{v}(ridx == 1, :)=[];
  end
  T = zeros(m,1);
  for i = 1:m
    Mi = cell(1,nvar);
    for v = 1:nvar
      Mi{v} = M{v}(1:m-1);
      M{v}(1:m-1)=[];
    end
    T(i,:) = feval(func, Mi{:});
  end
  Tori = mean(T, 1);
  Tori = Tori(ones(m,1), :);
  U = ((m-1) * (Tori-T));
  Var = (m-1) / m * sum((T-Tori).^2, 1);

  % Calculate standard error(s) of the functional parameter
  SE = sqrt(Var);

end
