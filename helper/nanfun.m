function Y = nanfun (func, X)

  % Private function file required for ibootci

  % Math functions, ignoring NaNs.
  [m,n] = size(X);
  if m > 1
    Y = zeros(1,n);
    for i = 1:n
      Y(i) = feval(func, X(~isnan(X(:,i)),i));
    end
  else
    Y = feval(func, X(~isnan(X)));
  end

end
