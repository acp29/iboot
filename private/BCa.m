function [m1, m2, SE, T] = BCa (nboot, func, data, bootstat, stat, alpha, S, paropt)

  % Private function file required for bootci

  % Redefine alpha as nominal coverage
  alpha = 1 - alpha;

  % Create distribution functions
  stdnormcdf = @(x) 0.5 * (1 + erf (x / sqrt(2)));
  stdnorminv = @(p) sqrt (2) * erfinv (2 * p-1);

  % Calculate bias correction z0
  z0 = stdnorminv (sum (bootstat < stat)/ nboot);

  % Use the Jackknife to calculate acceleration
  [SE, T, U] = jack (data, func, paropt);
  a = (1 / 6) * (sum(U.^3) / sum(U.^2)^(3/2));

  % Calculate BCa percentiles
  z1 = stdnorminv(0.5 * (1 + alpha));
  m1 = stdnormcdf(z0 + ((z0 + z1)/(1 - a * (z0 + z1))));
  z2 = stdnorminv(0.5 * (1 - alpha));
  m2 = stdnormcdf(z0 + ((z0 + z2)/(1 - a * (z0 + z2))));
  S.z0 = z0;
  S.a = a;

end
