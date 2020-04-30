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
  T = zeros(m,1);
  i = [1:m].';
  try
    pool = gcp('nocreate');
  catch
    pool = [];
  end
  if paropt.UseParallel && isoctave
    % Octave parallel computing
    parfun = @(i) parjack (x, func, nvar, m, i);
    T = pararrayfun(paropt.nproc, parfun, i, 'UniformOutput',true);
  elseif (paropt.UseParallel || ~isempty(pool)) && ~isoctave
    % Matlab parallel computing
    parfor i = 1:m
      M = cell(1,nvar);
      for v = 1:nvar
        M{v} = x{v};
        M{v}(i)=[];
      end
      T(i,:) = feval(func, M{:});
      M = [];
    end
  else
    % Octave or Matlab serial computing
    arrfun = @(i) parjack (x, func, nvar, m, i);
    T = arrayfun(arrfun, i, 'UniformOutput',true);
  end
  Tori = mean(T, 1);
  Tori = Tori(ones(m,1), :);
  U = ((m-1) * (Tori-T));
  Var = (m-1) / m * sum((T-Tori).^2, 1);

  % Calculate standard error(s) of the functional parameter
  SE = sqrt(Var);

end
