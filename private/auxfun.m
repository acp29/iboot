function T = auxfun (bootfun, nvar, varargin)

  % Private function file required for ibootci
  % Auxiliary function for block bootstrap
  X = varargin{1};
  Y = cat_blocks(nvar,X{:});
  T = bootfun(Y{:});

end
