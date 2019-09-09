%  Function File: iboottest2
%
%  Two-sample bootstrap test
%
%   p = iboottest2(nboot,{bootfun,x,y})
%   p = iboottest2(nboot,{bootfun,x,y},Name,Value)
%   [p,ci] = iboottest2(...)
%   [p,ci,S] = iboottest2(...)
%
%  Two sample bootstrap test for univariate data. The null hypothesis
%  is that the difference between the bootfun statistic calculated
%  for independent samples x and y is equal to zero. The test is
%  two-tailed.
%
%  See ibootci documentation for input argument definitions and
%  for Name-Value pairs. The following two input arguments differ
%  in their format:
%
%  'Weights': The 'Weights' option must be provided as a cell array:
%  the first cell should contain weights for x and the second cell
%  weights should contain weights for y. An empty cell signifies
%  that no weights will be used in the bootstrap for that sample.
%
%  'Strata': The same format as for weights.
%
%  'Clusters': The same format as for weights.
%
%  The syntax in this function code is known to be compatible with
%  recent versions of Octave (v3.2.4 on Debian 6 Linux 2.6.32) and
%  Matlab (v7.4.0 on Windows XP).
%
%  iboottest2 v1.4.0.0 (09/09/2019)
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


function [p,ci,S] = iboottest2(argin1,argin2,varargin)

  % Check and process iboottest2 input arguments
  nboot = argin1;
  bootfun = argin2{1};
  x = argin2{2};
  if all(size(x)>1)
    error('data must be a vector')
  end
  if size(x,2)>1
    x = x.';
  end
  y = argin2{3};
  if all(size(y)>1)
    error('data must be a vector')
  end
  if size(y,2)>1
    y = y.';
  end
  if numel(argin2)>3
    error('too many data variables provided')
  end
  iter = numel(nboot);
  if iter > 2
    error('Size of nboot exceeds maximum number of iterations supported by ibootci')
  end
  if iter==0
    B = 5000;
    C = 200;
    nboot = [B C];
  elseif iter==1
    B = nboot;
    C = 0;
    nboot = [B C];
  elseif iter==2
    B = nboot(1);
    C = nboot(2);
  end

  % Retireve some ibootci options
  options = varargin;
  alpha = 1+find(strcmpi('alpha',options));
  weights = 1+find(strcmpi('Weights',options));
  strata = 1+find(strcmpi('Strata',options));
  clusters = 1+find(strcmpi('Clusters',options));
  cellref = [];
  if ~isempty(alpha)
    try
      cellref = cat(2,cellref,[alpha-1,alpha]);
      alpha = options{alpha};
    catch
      alpha = 0.05;
      cellref(end-1:end)=[];
    end
  else
    alpha = 0.05;
  end
  if ~isa(alpha,'numeric') || numel(alpha)~=1
    error('The alpha value must be a numeric scalar value');
  end
  if (alpha <= 0) || (alpha >= 1)
    error('The alpha value must be a value between 0 and 1');
  end
  if ~isempty(weights)
    try
      cellref = cat(2,cellref,[weights-1,weights]);
      weights = options{weights};
    catch
      weights = {[],[]};
      cellref(end-1:end)=[];
    end
  else
    weights = {[],[]};
  end
  if ~isempty(strata)
    try
      cellref = cat(2,cellref,[strata-1,strata]);
      strata = options{strata};
    catch
      strata = {[],[]};
      cellref(end-1:end)=[];
    end
  else
    strata = {[],[]};
  end
  if ~isempty(clusters)
    try
      cellref = cat(2,cellref,[clusters-1,clusters]);
      clusters = options{clusters};
    catch
      clusters = {[],[]};
      cellref(end-1:end)=[];
    end
  else
    clusters = {[],[]};
  end
  options(cellref)=[];   % remove these evaluated options from the options array

  % Perform independent resampling from x and y
  state = warning;
  warning off;
  [~,bootstatX,SX] = ibootci(nboot,{bootfun,x},'Strata',strata{1},'Clusters',clusters{1},options{:});
  [~,bootstatY,SY] = ibootci(nboot,{bootfun,y},'Strata',strata{2},'Clusters',clusters{2},options{:});
  if C>0
    if ~isempty(weights{1})
      [~,bootstatX{1}] = ibootci(B,{bootfun,x},'alpha',SX.cal,'Weights',weights{1},'Strata',strata{1},options{:});
    end
    if ~isempty(weights{2})
      [~,bootstatY{1}] = ibootci(B,{bootfun,y},'alpha',SY.cal,'Weights',weights{2},'Strata',strata{2},options{:});
    end
  end
  warning(state);

  % Calculate differences between bootfun evaluated for bootstrap
  % sample sets derived from x and y
  if iscell(bootstatX)
    bootstatZ = cell(2,1);
    bootstatZ{1} = bootstatX{1} - bootstatY{1};
    bootstatZ{2} = bootstatX{2} - bootstatY{2};
  else
    bootstatZ = bootstatX - bootstatY;
  end

  % Create template settings structure and calculate the sample test statistic
  S = SX;
  T0 = SX.stat - SY.stat;
  S.stat = T0;           % assign correct sample test statistic in S
  S.alpha = alpha;       % reset alpha in S
  S.nboot = nboot;       % reset nboot in S
  S.n = SX.n + SY.n;     % calculate total sample size
  S.df = SX.df + SY.df;  % calculate degrees of freedom

  % Calculate confidence interval using ibootci
  [ci,bootstat,S,calcurve] = ibootci(bootstatZ, S);

  % If applicable, remove strata information from the output structure
  if ~isempty(strata{1})
    S = rmfield(S,{'SSb','SSw','ICC'});
  end
  S.strata = strata;
  S.clusters = clusters;
  S.weights = weights;

  % Calculate p-value using ibootp
  p = ibootp(0,bootstat,S,calcurve);

end
