function T = parjack (x, func, nvar, m, i)

  % Helper function file required for Octave Parallel Computing
  % capabilities of jack.

  M = cell(1,nvar);
  for v = 1:nvar
    M{v} = x{v};
    M{v}(i)=[];
  end
  T = feval(func, M{:});
  M = [];

end
