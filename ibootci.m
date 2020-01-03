%  Function File: ibootci
%
%  Two-sided bootstrap confidence intervals for small samples and  
%  samples with complex dependence structures
%
%  ci = ibootci(nboot,bootfun,...)
%  ci = ibootci(nboot,{bootfun,...},...,'alpha',alpha)
%  ci = ibootci(nboot,{bootfun,...},...,'type',type)
%  ci = ibootci(nboot,{bootfun,...},...,'type','stud','nbootstd',nbootstd)
%  ci = ibootci(nboot,{bootfun,...},...,'type','stud','stderr',stderr)
%  ci = ibootci(nboot,{bootfun,...},...,'weights',weights)
%  ci = ibootci(nboot,{bootfun,...},...,'strata',strata)
%  ci = ibootci(nboot,{bootfun,...},...,'cluster',clusters)
%  ci = ibootci(nboot,{bootfun,...},...,'block',blocksize)
%  ci = ibootci(nboot,{bootfun,...},...,'smooth',bandwidth)
%  ci = ibootci(nboot,{bootfun,...},...,'bootidx',bootidx)
%  ci = ibootci(nboot,{bootfun,...},...,'DEFF',state)
%  [ci,bootstat] = ibootci(...)
%  [ci,bootstat,S] = ibootci(...)
%  [ci,bootstat,S,calcurve] = ibootci(...)
%  [ci,bootstat,S,calcurve,bootidx] = ibootci(...)
%
%  ci = ibootci(nboot,bootfun,...) computes the 95% iterated (double)
%  bootstrap confidence interval of the statistic computed by bootfun.
%  nboot is a scalar, or vector of upto two positive integers indicating
%  the number of replicate samples for the first and second bootstraps.
%  bootfun is a function handle specified with @, or a string indicating
%  the function name. The third and later input arguments are data (column
%  vectors, or a matrix), that are used to create inputs for bootfun.
%  ibootci creates each first level bootstrap by block resampling from the
%  rows of the data argument(s) (which must be the same size) [1]. If a
%  positive integer for the number of second bootstrap replicates is
%  provided, then nominal central coverage of two-sided intervals is
%  calibrated to achieve second order accurate coverage by bootstrap
%  iteration and interpolation [2]. Linear interpolation of the empirical
%  cumulative distribution function of bootstat is then used to construct
%  two-sided confidence intervals [3]. The resampling method used 
%  throughout is balanced resampling [4]. Default values for the number of 
%  first and second bootstrap replicate sample sets in nboot are 5000 and 
%  200 respectively. Note that this calibration procedure does not apply 
%  to Studentized-type intervals or Cluster or Weights bootstrap options.
%
%  ci = ibootci(nboot,{bootfun,...},...,'alpha',alpha) computes the
%  iterated bootstrap confidence interval of the statistic defined by the
%  function bootfun with coverage 100*(1-alpha)%, where alpha is a scalar
%  value between 0 and 1. bootfun and the data that ibootci passes to it
%  are contained in a single cell array. The default value of alpha is
%  0.05 corresponding to intervals with a coverage of 95% confidence.
%
%  ci = ibootci(nboot,{bootfun,...},...,'type',type) computes the bootstrap
%  confidence interval of the statistic defined by the function bootfun.
%  type is the confidence interval type, chosen from among the following:
%    'per' or 'percentile': Percentile method (Default). 
%    'stud' or 'student': Studentized (bootstrap-t) confidence interval.
%    The bootstrap-t method includes an additive correction to stabilize
%    the variance when the sample size is small [6].
%    'cper': Bias-corrected percentile method [7].
%
%  ci = ibootci(nboot,{bootfun,...},...,'type','stud','nbootstd',nbootstd)
%  computes the Studentized bootstrap confidence interval of the statistic
%  defined by the function bootfun. The standard error of the bootstrap
%  statistics is estimated using bootstrap, with nbootstd bootstrap data
%  samples. nbootstd is a positive integer value. The default value of
%  nbootstd is 200. Setting nbootstd overides the second element in nboot.
%  The nbootstd argument is ignored when the interval type is set to
%  anything other than Studentized (bootstrap-t) intervals.
%
%  ci = ibootci(nboot,{bootfun,...},...,'type','stud','stderr',stderr) 
%  computes the studentized bootstrap confidence interval of statistics 
%  defined by the function bootfun. The standard error of the bootstrap 
%  statistics is evaluated by the function stderr. stderr is a function 
%  handle. stderr takes the same arguments as bootfun and returns the 
%  standard error of the statistic computed by bootfun.
%
%  ci = ibootci(nboot,{bootfun,...},...,'weights',weights) specifies
%  observation weights. weights must be a vector of non-negative numbers.
%  The dimensions of weights must be equal to that of the non-scalar input
%  arguments to bootfun. The weights are used as bootstrap sampling
%  probabilities. Note that weights are not implemented for Studentized-
%  type intervals or bootstrap iteration.
%
%  ci = ibootci(nboot,{bootfun,...},...,'strata',strata) specifies a
%  vector containing numeric identifiers of strata. The dimensions of
%  strata must be equal to that of the non-scalar input arguments to
%  bootfun. Bootstrap resampling is stratified so that every stratum is
%  represented in each bootstrap test statistic [5]. If weights are also
%  provided then they are within-stratum weights; the weighting of
%  individual strata depends on their respective sample size.
%
%  ci = ibootci(nboot,{bootfun,...},...,'cluster',clusters) specifies
%  a vector containing numeric identifiers for clusters. Whereas strata
%  are fixed, clusters are resampled. This is achieved by two-stage
%  bootstrap resampling of residuals with shrinkage correction [5,8,9].
%  If a matrix is provided defining additional levels of subsampling in
%  a hierarchical data model, then level two cluster means are computed
%  and resampled. This option is not compatible with bootstrap iteration.
%
%  ci = ibootci(nboot,{bootfun,...},...,'block',blocksize) specifies
%  a positive integer defining the block length for block bootstrapping
%  data with serial dependence (e.g. stationary time series). The
%  algorithm uses circular, overlapping blocks. Intervals are constructed
%  without standardization making them equivariant under monotone
%  transformations [10]. The double bootstrap resampling and calibration
%  procedure makes interval coverage less sensitive to block length [11].
%  If the blocksize is set to 'auto' (recommended), the block length is
%  calculated automatically. Note that balanced resampling is not
%  maintained for block bootstrap. Block bootstrap can also be used for
%  regression of time series data by combining it with pairs bootstrap
%  (i.e. by providing x and y vectors as data variables).
%
%  ci = ibootci(nboot,{bootfun,...},...,'smooth',bandwidth) applies
%  additive random Gaussian noise of the specified bandwidth to the 
%  bootstrap sample sets before evaluating bootfun. Recommended usage 
%  for univariate data, is to run bootstrap using the appropriate 
%  resampling method first without smoothing, and return the output 
%  structure S. Then re-run ibootci and apply smoothing with the 
%  bandwidth set to the standard error (S.SE), which will be of the 
%  order n^(-1/2), suitable for improving coverage of confidence
%  intervals when sample sizes are small [12]. For d-dimensional 
%  multivariate data, the bandwidth can be an 1-by-d vector of 
%  bandwidths, or a d-by-d covariance matrix. 
%
%  ci = ibootci(nboot,{bootfun,...},...,'bootidx',bootidx) performs
%  bootstrap computations using the indices from bootidx for the first
%  bootstrap.
%
%  ci = ibootci(nboot,{bootfun,...},...,'DEFF',state) estimates the 
%  design effect by resampling. State can be 'on' or 'off'. Default is
%  'off'.
%
%  [ci,bootstat] = ibootci(...) also returns the bootstrapped statistic
%  computed for each of the bootstrap replicate samples sets. If only
%  a single bootstrap is requested, bootstat will return a vector: each
%  column of bootstat contains the result of applying bootfun to one
%  replicate sample from the first bootstrap. If a function handle is 
%  passed to stderr for Studentised bootstrap intervals, bootstat will 
%  return a cell array containing the statistics computed using bootfun 
%  and stderr in the first and second cells respectively. If bootstrap 
%  iteration is requested, bootstat will return a cell array containing 
%  the statistics computed by bootfun in the first and second bootstrap.
%  For the second boostrap, each column of bootstat contains the
%  results of applying bootfun to each replicate sample from the second
%  bootstrap for one replicate sample from the first bootstrap.
%
%  [ci,bootstat,S] = ibootci(...) also returns a structure containing
%  the settings used in the bootstrap and the resulting statistics
%  including the (double) bootstrap bias and standard error.
%
%  The output structure S contains the following fields:
%    bootfun: Function name or handle used to calculate the test statistic
%    nboot: The number of first (and second) bootstrap replicate samples
%    nvar: Number of data variables
%    n: The length of data variable(s) (and no. of clusters if applicable)
%    type: Type of confidence interval (per or stud)
%    alpha: Desired alpha level
%    coverage: Central coverage of the confidence interval
%    cal: Nominal alpha level from calibration
%    z0: Bias correction (0 if type is not 'cper')
%    bandwidth: Bandwidth for smooth bootstrap (Gaussian/Student-t kernel)
%    xcorr: Autocorrelation coefficients (maximum 99 lags)
%    ICC: Intraclass correlation coefficient - one-way random, ICC(1,1)
%    DEFF: Design effect estimated by resampling (if requested)
%    stat: Sample test statistic calculated by bootfun
%    bias: Bias of the test statistic
%    bc_stat: Bias-corrected test statistic
%    SE: Bootstrap standard error of the test statistic
%    ci: Bootstrap confidence interval of the test statistic
%    prct: Percentiles used to generate confidence intervals (proportion)
%    weights: Argument supplied to 'Weights' (empty if none provided)
%    strata: Argument supplied to 'Strata' (empty if none provided)
%    clusters: Argument supplied to 'Clusters' (empty if none provided)
%    blocksize: Length of overlapping blocks (empty if none provided)
%
%  [ci,bootstat,S,calcurve] = ibootci(...) also returns the calibration
%  curve for central coverage. The first column is nominal coverage and
%  the second column is actual coverage.
%
%  [ci,bootstat,S,calcurve,bootidx] = ibootci(...) also returns bootidx,
%  a matrix of indices from the first bootstrap.
%
%  Computations of confidence intervals can be accelerated by starting 
%  a parallel pool prior to executing ibootci. This is particularly useful
%  for the calculation of p-values using ibootci wrapper functions.
%  Parallel usage is supported in MATLAB for most features except for
%  balanced resampling during the first bootstrap.
%
%  Bibliography:
%  [1] Efron, and Tibshirani (1993) An Introduction to the
%        Bootstrap. New York, NY: Chapman & Hall
%  [2] Hall, Lee and Young (2000) Importance of interpolation when
%        constructing double-bootstrap confidence intervals. Journal
%        of the Royal Statistical Society. Series B. 62(3): 479-491
%  [3] Efron (1981) Censored data and the bootstrap. JASA
%        76(374): 312-319
%  [4] Davison et al. (1986) Efficient Bootstrap Simulation.
%        Biometrika, 73: 555-66
%  [5] Davison and Hinkley (1997) Bootstrap Methods and their
%        application. Chapter 3: pg 97-100
%  [6] Polansky (2000) Stabilizing bootstrap-t confidence intervals
%        for small samples. Can J Stat. 28(3):501-516
%  [7] Efron (1982) The jackknife, the bootstrap and other resampling
%        plans. CBMS 38, SIAM-NSF.
%  [8] Gomes et al. (2012) Developing appropriate methods for cost-
%        effectiveness analysis of cluster randomized trials.
%        Medical Decision Making. 32(2): 350-361
%  [9] Ng, Grieve and Carpenter (2013) Two-stage nonparametric
%        bootstrap sampling with shrinkage correction for clustered
%        data. The Stata Journal. 13(1): 141-164
%  [10] Gotze and Kunsch (1996) Second-Order Correctness of the Blockwise
%        Bootstrap for Stationary Observations. The Annals of Statistics.
%        24(5):1914-1933
%  [11] Lee and Lai (2009) Double block bootstrap confidence intervals
%        for dependent data. Biometrika. 96(2):427-443
%  [12] Polansky and Schucany (1997) Kernel Smoothing to Improve Bootstrap 
%        Confidence Intervals. J R Statist Soc B. 59(4):821-838 
%
%  Example 1: Two alternatives for 95% confidence intervals for the mean
%    >> y = randn(20,1);
%    >> ci = ibootci([5000 200],@mean,y);
%    >> ci = ibootci([5000 200],{@mean,y},'alpha',0.05);
%
%  Example 2: 95% confidence intervals for the means of paired/matched data
%    >> y1 = randn(20,1);
%    >> y2 = randn(20,1);
%    >> [ci1,bootstat,S,calcurve,bootidx] = ibootci([5000 200],{@mean,y1});
%    >> ci2 = ibootci([5000 200],{@mean,y2},'bootidx',bootidx);
%
%  Example 3: 95% confidence intervals for the correlation coefficient
%    >> z = mvnrnd([2,3],[1,1.5;1.5,3],20);
%    >> x = z(:,1); y = z(:,2);
%    >> corrcoef = @(X,Y) diag(corr(X,Y)).';
%    >> ci = ibootci([5000 200],{corrcoef,x,y});
%  Note that this is much faster than:
%    >> ci = ibootci([5000 200],{@corr,x,y});
%
%  Example 4: 95% confidence interval for the weighted arithmetic mean
%    >> y = randn(20,1);
%    >> w = [ones(5,1)*10/(20*5);ones(15,1)*10/(20*15)];
%    >> [ci,bootstat,S] = ibootci([5000,200],{'mean',y},'alpha',0.05);
%    >> ci = ibootci(5000,{'mean',y},'alpha',S.cal,'Weights',w);
%
%  Example 5: 95% confidence interval for the median by smoothed bootstrap
%  (requires the smoothmedian function available at Matlab Central File 
%  Exchange)
%    >> y = randn(20,1);
%    >> ci = ibootci([5000 200],@smoothmedian,y);
%
%  Example 6: 95% confidence interval for the 25% trimmed (or 
%  interquartile) mean
%    >> y = randn(20,1);
%    >> func = @(x) trimmean(x,50)
%    >> ci = ibootci([5000 200],func,y);
%
%  The syntax in this function code is known to be compatible with
%  recent versions of Octave (v3.2.4 on Debian 6 Linux 2.6.32) and
%  Matlab (v7.4.0 on Windows XP).
%
%  ibootci v2.8.4.2 (02/01/2020)
%  Author: Andrew Charles Penn
%  https://www.researchgate.net/profile/Andrew_Penn/
%
%  Cite as:
%  Andrew Penn (2019). ibootci (https://www.github.com/acp29/iboot), GitHub.
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


function [ci,bootstat,S,calcurve,idx] = ibootci(argin1,argin2,varargin)

  % Evaluate the number of function arguments
  if nargin<2
    error('Too few input arguments');
  end
  if nargout>5
   error('Too many output arguments');
  end

  % Assign input arguments to function variables
  if isstruct(argin2)
    %  Usage (undocumented):
    %
    %  ci = ibootci(bootstat,S) produces (calibrated) confidence intervals
    %  for the bootstrap replicate sample set statistics provided in bootstat.
    %  This usage also requires a complete settings structure.
    %
    bootstat = argin1;
    if iscell(bootstat)
      T1 = bootstat{1};
      T2 = bootstat{2};
      B = numel(T1);
      if size(T2,1) > 1
        C = size(T2,1);
      else
        C = 0;
        stderr = [];
      end
    else
      T1 = bootstat;
      B = numel(T1);
      C = 0;
    end
    S = argin2;
    nboot = [B,C];
    if B ~= S.nboot(1)
      error('the dimensions of bootstat are inconsistent with S.nboot')
    end
    bootfun = S.bootfun;
    data = [];
    ori_data = data;
    n = S.n(1);
    nvar = S.nvar;
    idx = [];
    T0 = S.stat;
    bandwidth = S.bandwidth;
    weights = S.weights;
    strata = S.strata;
    clusters = S.clusters;
    blocksize = S.blocksize;
    deff = 'off';
    type = S.type;
    S.coverage = 1-S.alpha;
    alpha = S.coverage;       % convert alpha to coverage
    opt = struct;
    opt.weights = weights;
    opt.strata = strata;
    opt.clusters = clusters;
    opt.blocksize = blocksize;
    opt.bandwidth = bandwidth;
    % Perform calibration (if applicable)
    if C>0 && any(strcmpi(type,{'per','percentile','cper'}))
      U = zeros(1,B);
      for h = 1:B
        U(h) = interp_boot2(T2(:,h),T0,C);
      end
    end
    
  elseif ~iscell(argin2)
    % Normal usage without options
    nboot = argin1;
    bootfun = argin2;
    data = varargin;
    ori_data = data;
    alpha = 0.05;
    idx = [];
    weights = [];
    strata = [];
    clusters = [];
    blocksize = [];
    bandwidth = [];
    nbootstd = [];
    deff = 'off';
    type = 'per';
    T1 = [];  % Initialize bootstat variable

  else
    % Normal usage with options
    % Evaluate option input arguments
    nboot = argin1;
    bootfun = argin2{1};
    data = {argin2{2:end}};
    options = varargin;
    alpha = 1+find(strcmpi('alpha',options));
    type = 1+find(strcmpi('type',options));
    nbootstd = 1+find(strcmpi('nbootstd',options));
    stderr = 1+find(strcmpi('stderr',options));
    weights = 1+find(strcmpi('weights',options));
    strata = 1+find(cellfun(@(options) any(strcmpi({'strata','stratum','stratified'},options)),options));
    clusters = 1+find(cellfun(@(options) any(strcmpi({'clusters','cluster'},options)),options));
    blocksize = 1+find(cellfun(@(options) any(strcmpi({'block','blocks','blocksize'},options)),options));
    bandwidth = 1+find(cellfun(@(options) any(strcmpi({'smooth','smoothing','bandwidth'},options)),options));
    bootidx = 1+find(strcmpi('bootidx',options));
    deff = 1+find(strcmpi('deff',options));
    if ~isempty(alpha)
      try
        alpha = options{alpha};
      catch
        alpha = 0.05;
      end
    else
      alpha = 0.05;
    end
    if ~isempty(type)
      try
        type = options{type};
      catch
        type = 'per';
      end
    else
      type = 'per';
    end
    if any(strcmpi(type,{'stud','student'}))
      if ~isempty(nbootstd)
        try
          nbootstd = options{nbootstd};
          nboot = [nboot(1) nbootstd];
        catch
          nbootstd = 200;
          nboot = [nboot(1),nbootstd];
        end
      else
        nbootstd = [];
      end
    else
      nbootstd = [];
    end
    if ~isempty(stderr)
      try
        stderr = options{stderr};
      catch
        stderr = [];
      end
    else
      stderr = [];
    end
    if ~isempty(weights)
      try
        weights = options{weights};
      catch
        weights = [];
      end
    else
      weights = [];
    end
    if ~isempty(strata)
      try
        strata = options{strata};
      catch
        strata = [];
      end
    else
      strata = [];
    end
    if ~isempty(clusters)
      try
        clusters = options{clusters};
        if ~isempty(clusters)
          strata = clusters;
        end
      catch
        clusters = [];
      end
    else
      clusters = [];
    end
    if ~isempty(blocksize)
      try
        blocksize = options{blocksize};
      catch
        blocksize = [];
      end
    else
      blocksize = [];
    end
    if ~isempty(bandwidth)
      try
        bandwidth = options{bandwidth};
        if ~all(bandwidth(:))
          bandwidth = [];
        end
      catch
        bandwidth = [];
      end
    else
      bandwidth = [];
    end
    if ~isempty(bootidx)
      try
        idx = options{bootidx};
      catch
        error('Could not find bootidx')
      end
      if size(data{1},1) ~= size(idx,1)
        error('Dimensions of data and bootidx are inconsistent')
      end
      % Set nboot(1) according to the size of bootidx
      nboot(1) = size(idx,2);
    else
      idx = [];
    end
    if ~isempty(deff)
      try
        deff = options{deff};
      catch
        deff = 'off';
      end
    else
      deff = 'off';
    end
    T1 = [];  % Initialize bootstat variable
  end

  if isempty(T1)
    % Evaluate function variables
    iter = numel(nboot);
    if iter > 2
      error('Size of nboot exceeds maximum number of iterations supported by ibootci')
    end
    if ~isa(nboot,'numeric')
      error('nboot must be numeric');
    end
    if any(nboot~=abs(fix(nboot)))
      error('nboot must contain positive integers')
    end
    if ~isa(alpha,'numeric') || numel(alpha)~=1
      error('The alpha value must be a numeric scalar value');
    end
    if (alpha <= 0) || (alpha >= 1)
      error('The alpha value must be a value between 0 and 1');
    end
    if ~any(strcmpi(type,{'per','percentile','stud','student','cper'}))
      error('The type of bootstrap must be either per or stud');
    end
    if ~isa(alpha,'numeric')
      error('The smoothing bandwidth(s) must be numeric')
    end
    if sum(strcmpi(deff,{'on','off'})) < 1
      error('The deff input argument must be set to ''on'' or ''off''')
    end
    if ~isempty(bandwidth)
      if strcmpi(deff,'on')
        error('Design effect calculation not compatible with smoothing')
      end
    end
    
    % Evaluate data input
    nvar = size(data,2);
    if (min(size(data{1}))>1)
      if (nvar == 1)
        nvar = size(data{1},2);
        data = num2cell(data{1},1);
        matflag = 1;   % Flag for matrix input set to True
        runmode = 'slow';
      else
        error('Multiple data input arguments must be provided as vectors')
      end
    else
      matflag = 0;
      runmode = [];
    end
    varclass = zeros(1,nvar);
    rows = zeros(1,nvar);
    cols = zeros(1,nvar);
    for v = 1:nvar
      varclass(v) = isa(data{v},'double');
      if all(size(data{v})>1) && (v > 1)
        error('Vector input arguments must be the same size')
      end
      rows(v) = size(data{v},1);
      cols(v) = size(data{v},2);
    end
    if ~all(varclass)
      error('Data variables must be double precision')
    end
    if any(rows~=rows(1)) || any(cols~=cols(1))
      error('The dimensions of the data are not consistent');
    end
    rows = rows(1);
    cols = cols(1);
    if max(rows,cols) == 1
      error('Cannot bootstrap scalar values');
    elseif cols>1
      % Transpose row vector data
      n = cols;
      for v = 1:nvar
        data{v} = data{v}.';
      end
    else
      n = rows;
    end
    ori_data = data; % Make a copy of the data
    if ~isempty(strata) || ~isempty(clusters)
      if size(clusters,1) ~= size(data{1},1)
        error('Dimensions of clusters are inconsistent with the data')
      end
      if size(strata,1) ~= size(data{1},1)
        error('Dimensions of strata are inconsistent with the data')
      end
      while size(strata,2) > 1
        % Calculate strata means for resampling more than two nested levels
        % Picquelle and Mier (2011) Fisheries Research 107(1-3):1-13
        [data,strata] = unitmeans(data,strata,nvar);
      end
      if ~isempty(clusters)
        clusters = strata;
      end
      if numel(unique(strata)) == 1
        strata = []; % Cannot perform stratified resampling
      else
        n = size(strata,1);
      end
      % Sort strata and data vectors so that strata components are grouped
      [strata,I] = sort(strata);
      for v = 1:nvar
        data{v} = data{v}(I);
      end
      ori_data = data; % recreate copy of the data
    end
    if isempty(weights)
      weights = ones(n,1);
    else
      if size(weights,2)>1
        % Transpose row vector weights
        weights = weights.';
      end
      if ~all(size(weights) == [n,1])
        error('The weights vector is not the same dimensions as the data');
      end
    end
    if any(weights<0)
      error('Weights must be a vector of non-negative numbers')
    end
    % Evaluate bootfun
    if ischar(bootfun)
      % Convert character string of a function name to a function handle
      bootfun = str2func(bootfun);
    end
    if ~isa(bootfun,'function_handle')
      error('bootfun must be a function name or function handle');
    end
    try
      if matflag > 0
        temp = list2mat(data{:});
        T0 = feval(bootfun,temp);
      else
        T0 = feval(bootfun,data{:});
      end
    catch
      error('An error occurred while trying to evaluate bootfun with the input data');
    end
    if isinf(T0) || isnan(T0)
      error('bootfun returns a NaN or Inf')
    end
    if max(size(T0))>1
      error('Column vector inputs to bootfun must return a scalar');
    end
    % Minimal simulation to evaluate bootfun with matrix input arguments
    M = cell(1,nvar);
    for v = 1:nvar
      x = data{v};
      if v == 1
        simidx = ceil(n.*rand(n,2));  % For compatibility with R2007
      end
      M{v} = x(simidx);
    end
    if isempty(runmode)
      try
        sim = feval(bootfun,M{:});
        if size(sim,1)>1
          error('Invoke catch statement');
        end
        runmode = 'fast';
      catch
        warning('ibootci:slowMode',...
                'Slow mode. Faster if matrix input arguments to bootfun return a row vector.')
        runmode = 'slow';
      end
      if ~isempty(stderr)
        if ischar(bootfun)
          % Convert character string of a function name to a function handle
          stderr = str2func(stderr);
        end
        if ~isa(stderr,'function_handle')
          error('stderr must be a function name or function handle');
        end
        if strcmpi(runmode,'fast')
          try
            simSE = feval(stderr,M{:});
            if size(simSE,1)>1
              error('Invoke catch statement');
            end
              runmode = 'fast';
          catch
            warning('ibootci:slowMode',...
                    'Slow mode. Faster if matrix input arguments to stderr return a row vector.')
            runmode = 'slow';
          end
        end
      end
    end

    % Set the bootstrap sample sizes
    if iter==0
      B = 5000;
      C = 200;
      nboot = [B,C];
    elseif iter==1
      B = nboot;
      C = 0;
      nboot = [B,C];
    elseif iter==2
      B = nboot(1);
      C = nboot(2);
    end
    if ~isempty(stderr)
      if C > 0
        error('bootstrap iteration not compatible with the stderr option')
      end
      nboot(2) = 0;
      C = nboot(2);
    end
    if isempty(nbootstd) && isempty(stderr) && (C==0) && any(strcmpi(type,{'stud','student'})) 
      error('Studentized (bootstrap-t) intervals require bootstrap interation or stderr')
    end
    if C>0 && ~any(strcmpi(type,{'stud','student'}))
      if (1/min(alpha,1-alpha)) > (0.5*C)
        error('ibootci:extremeAlpha',...
             ['The calibrated alpha is too extreme for calibration so the result will be unreliable. \n',...
              'Try increasing the number of replicate samples in the second bootstrap.\n',...
              'If the problem persists, the original sample size may be inadequate.\n']);
      end
      if any(diff(weights))
        error('Weights are not implemented for iterated bootstrap.');
      end
    end

    % Initialize output structure
    S = struct;
    S.bootfun = bootfun;
    S.nboot = nboot;
    S.nvar = nvar;
    S.n = n;
    S.type = type;
    S.alpha = alpha;
    S.coverage = 1-alpha;
    S.bandwidth = bandwidth;

    % Update bootfun for matrix data input argument
    if matflag > 0
      bootfun = @(varargin) bootfun(list2mat(varargin{:}));
    end
    
    % Prepare bandwidth variable and correlation matrix
    if ~isempty(bandwidth)
      if C>0 
        error('Smoothing not available for iterated bootstrap')
      end
      if (min(size(bandwidth)) > 1)
        % Do nothing, bandwidth is already a covariance matrix
      else
        % Format bandwidth input (scalar or vector) as covariance matrix
        bandwidth = diag(bandwidth.^2);
      end
      % Extract bandwidth and correlation matrix from covariance matrix
      D = diag(sqrt(diag(bandwidth))).'; % standard deviations
      R = inv(D) * bandwidth * inv(D);   % correlation matrix 
      bandwidth = diag(D).';             % set standard deviation as bandwidth 
      % Calculate degrees of freedom for Student-t distribution
      if numel(S.n)>1
        df = S.n(1)-S.n(2);
      else
        df = S.n(1) - 1;
      end
    else 
      R = [];
      df = [];
    end

    % Convert alpha to coverage level (decimal format)
    alpha = 1-alpha;

    % Prepare for cluster resampling (if applicable)
    if ~isempty(clusters)
      if ~isempty(blocksize)
        error('Incompatible combination of options.')
      end
      if any(strcmpi(type,{'stud','student'}))
        error('Bootstrapping clustered data is not implemented with bootstrap-t intervals.')
      end
      if C > 0
        error('Bootstrapping clustered data is not implemented with bootstrap iteration.')
      end
      if nargout > 4
        error('No bootidx for two-stage resampling of clustered data.')
      end
      % Redefine data as intracluster residuals
      % Residuals will undergo balanced resampling with replacement
      % Ordinary resampling with replacement is used for cluster means
      [mu,Z,K,g] = clustmean(data,clusters,nvar);
      data = Z;
      S.n(2) = K; % S.n is [number of observations, number of clusters]
      bootfun = @(varargin) bootclust(bootfun,K,g,runmode,mu,varargin);
    end

    % Prepare for block resampling (if applicable)
    if ~isempty(blocksize)
      if ~isempty(clusters) || ~isempty(strata) || any(diff(weights))
        error('Incompatible combination of options')
      end
      if strcmpi(blocksize,'auto')
        blocksize = round(n^(1/3));  % set block length to ~ n^(1/3)
      end
      data = split_blocks(data,blocksize);
      bootfun = @(varargin) auxfun(bootfun,S.nvar,varargin);
      nvar = S.nvar * blocksize;
    end

    % Prepare bootstrap options stucture
    opt = struct;
    opt.weights = weights;
    opt.strata = strata;
    opt.clusters = clusters;
    opt.blocksize = blocksize;
    opt.bandwidth = bandwidth;
    opt.R = R;
    opt.df = df;
    opt.stderr = stderr;
    opt.runmode = runmode;
    
    % Perform bootstrap
    % Bootstrap resampling
    if isempty(idx)
      if (nargout < 5)
        [T1, T2, U] = boot1 (data, nboot, n, nvar, bootfun, T0, S, opt);
      else
        [T1, T2, U, idx] = boot1 (data, nboot, n, nvar, bootfun, T0, S, opt);
      end
      if any(strcmpi(type,{'stud','student'})) && ~isempty(stderr) 
        % When a function is provided in stderr, the U variable contains 
        % standard errors of the bootstrap samples 
        SE1 = U;
      end
    else
      X1 = cell(1,nvar);
      for v = 1:nvar
        X1{v} = data{v}(idx);
      end
      if C>0
        if ~isempty(strata)
          [SSb, SSw, K, g] = sse_calc (data, strata, nvar);
        else
          g = ones(n,1);
        end
        T2 = zeros(C,B);
        U = zeros(1,B);
        for h = 1:B
          x1 = cell(1,nvar);
          for v = 1:nvar
            x1{v} = X1{v}(:,h);
          end
          [U(h), T2(:,h)] = boot2 (x1, nboot, n, nvar, bootfun, T0, g, S, opt);
        end
      else
        T2 = [];
      end
      if ~isempty(bandwidth)
        % Apply smoothing using Gaussian/Student-t kernel
        noise = bsxfun(@times,mvtrnd(R,df,n),bandwidth);
        for v = 1:nvar
          X1{v} = X1{v} + noise(:,v);
        end
      end
      switch lower(runmode)
        case {'fast'}
          T1 = feval(bootfun,X1{:});
          if any(strcmpi(type,{'stud','student'})) && ~isempty(stderr)
            SE1 = feval(stderr,X1{:});
          end
        case {'slow'}
          T1 = zeros(1,nboot(1));
          for h = 1:nboot(1)
            x1 = cellfun(@(X1)X1(:,h),X1,'UniformOutput',false);
            T1(h) = feval(bootfun,x1{:});
            if any(strcmpi(type,{'stud','student'})) && ~isempty(stderr)
              SE1(:,h) = feval(stderr,x1{:});
            end
          end
      end
    end
    % Assign data to bootstat
    if isempty(T2)
      if any(strcmpi(type,{'stud','student'})) && ~isempty(stderr) 
        bootstat = cell(2,1);
        bootstat{1} = T1;
      else 
        bootstat = T1;
      end
    else
      bootstat = cell(2,1);
      bootstat{1} = T1;
      bootstat{2} = T2;
    end
  end

  % Calculate statistics for the first bootstrap sample set
  if C>0
    % Double bootstrap bias estimation
    % See Davison and Hinkley (1997) pg 103-107
    b = mean(T1)-T0;
    c = mean(mean(T2))-2*mean(T1)+T0;
    bias = b-c;
    % Double bootstrap multiplicative correction of the variance
    SE = sqrt(var(T1,0)^2 / mean(var(T2,0)));
  else
    % Single bootstrap bias estimation
    bias = mean(T1)-T0;
    % Single bootstrap variance estimation
    SE = std(T1,0);
  end

  % Calibrate central two-sided coverage
  if C>0 && any(strcmpi(type,{'per','percentile','cper'}))
    % Create a calibration curve
    V = abs(2*U-1);
    [calcurve(:,2),calcurve(:,1)] = empcdf(V,1);
    alpha = interp1(calcurve(:,2),calcurve(:,1),alpha);
  else
    calcurve = [];
  end
  S.cal = 1-alpha;

  % Check the nominal central coverage
  if (S.cal == 0)
    warning('ibootci:calibrationHitEnd',...
            ['The calibration of alpha has hit the ends of the bootstrap distribution \n',...
             'and may be unreliable. Try increasing the number of replicate samples for the second \n',...
             'bootstrap. If the problem persists, the original sample size may be inadequate.\n']);
  end

  % Construct confidence interval (with calibrated central coverage)
  switch lower(type)
    case {'per','percentile'}
      % Percentile
      m1 = 0.5*(1+alpha);
      m2 = 0.5*(1-alpha);
      S.z0 = 0;
    case 'cper'
      % Bias corrected intervals
      [m1,m2,S] = BC(B,bootfun,data,T1,T0,alpha,S);
    case {'stud','student'}
      % Bootstrap-t
      m1 = 0.5*(1-alpha);
      m2 = 0.5*(1+alpha);
      S.z0 = 0;
  end

  % Linear interpolation for interval construction
  if any(strcmpi(type,{'stud','student'}))

    % Use bootstrap-t method with variance stabilization for small samples
    % Polansky (2000) Can J Stat. 28(3):501-516
    if C>0
      se = std(T1,0);
      SE1 = std(T2,0);
    else
      if ~isempty(stderr)
        % Use stderr function if provided
        se = feval(stderr,data{:});
      else
        % The standard errors should already be defined
        se = S.SE;
        SE1 = bootstat{2};
      end
      SE = se;
      % Put standard errors for into second cell of bootstat output 
      bootstat{2} = SE1;
    end
    if isempty(bandwidth) 
      a = n^(-3/2) * se;
    else
      a = 0;
    end

    % Calculate Studentized statistics
    T = (T1-T0)./(SE1+a);
    [cdf,T] = empcdf(T,1);

    % Calculate intervals from empirical distribution of the Studentized bootstrap statistics
    UL = T0 - se * interp1(cdf,T,m1,'linear','extrap');
    LL = T0 - se * interp1(cdf,T,m2,'linear','extrap');
    ci = [LL;UL];

  else

    % Calculate interval for percentile method
    [cdf,t1] = empcdf(T1,1);
    UL = interp1(cdf,t1,m1,'linear','extrap');
    LL = interp1(cdf,t1,m2,'linear','extrap');
    ci = [LL;UL];

  end

  % Check the confidence interval limits
  if (m2 < cdf(2)) || (m1 > cdf(end-1))
    warning('ibootci:intervalHitEnd',...
            ['The confidence interval has hit the end(s) of the bootstrap distribution \n',...
             'and may be unreliable. Try increasing the number of replicate samples in the second \n',...
             'bootstrap. If the problem persists, the original sample size may be inadequate.\n']);
  end

  % Analysis of dependence structure of the data
  % Also re-sort data to match original input data
  if (nargout>2)
    if ~isempty(data)
        
      % Examine dependence structure of each variable by autocorrelation
      if ~isempty(ori_data)
        S.xcorr = zeros(2*min(S.n(1),99)+1,S.nvar);
        for v = 1:S.nvar
          S.xcorr(:,v) = xcorr(ori_data{v},min(S.n(1),99),'coeff');
        end
        S.xcorr(1:min(S.n(1),99),:) = [];
      end  
        
      % Calculate intraclass correlation coefficient (ICC) for each variable
      %  - Smeeth and Ng (2002) Control Clin Trials. 23(4):409-21
      %  - Huang (2018) Educ Psychol Meas. 78(2):297-318
      %  - McGraw & Wong (1996) Psychological Methods. 1(1):30-46
      if ~isempty(strata) || ~isempty(clusters)
        % Intraclass correlation coefficient (based on the mean)
        % One-way random, single measures ICC(1,1)
        if ~isempty(strata)
          groups = strata;
        elseif ~isempty(clusters)
          groups = clusters;
        end
        [SSb, SSw, K, g, MSb, MSw, dk] = sse_calc (ori_data, groups, nvar);
        S.ICC = (MSb-MSw)/(MSb+(dk-1)*MSw);
      else
        S.ICC = NaN;
      end

      % Estimate the design effect by resampling
      % Ratio of variance to that calculated by simple random sampling (SRS)
      if matflag > 0
        bootfun = @(varargin) S.bootfun(list2mat(varargin{:}));
      else 
        bootfun = S.bootfun;
      end
      opt.weights = ones(n,1);
      opt.strata = [];
      opt.clusters = [];
      opt.blocksize = [];
      opt.bandwidth = [];
      opt.R = [];
      opt.df = [];
      if strcmpi(deff,'on')
        if ~isempty(clusters)
          [SRS1,SRS2] = boot1(ori_data,[B,min(B,200)],S.n(1),S.nvar,bootfun,T0,S,opt);
        else
          [SRS1,SRS2] = boot1(ori_data,S.nboot,S.n(1),S.nvar,bootfun,T0,S,opt);
        end
        if (C > 0) || ~isempty(clusters)
          SRSV = var(SRS1,0)^2 / mean(var(SRS2,0));
        else
          SRSV = var(SRS1,0);
        end
        S.DEFF = SE^2/SRSV;
      else
        S.DEFF = [];
      end

      % Re-sort variables to match input data
      if ~isempty(strata) || ~isempty(clusters)
        [sorted,J] = sort(I);
        strata = strata(J,:);
        if ~isempty(idx)
          idx = idx(J,:);
        end
      end
      if ~isempty(blocksize) && ~isempty(idx)
        temp = cell(n,1);
        for i = 1:n
          temp{i} = bsxfun(@plus,(0:blocksize-1)'.*ones(1,B),...
                                 ones(blocksize,1).*idx(i,:));
        end
        idx = cell2mat(temp);
        idx(n+1:end,:) = [];
        idx(idx>n) = idx(idx>n)-n;
      end

    else
      S.ICC = [];
      S.DEFF = [];
      S.xcorr = [];
    end

    % Complete output structure
    S.stat = T0;         % Sample test statistic
    S.bias = bias;       % Bias of the test statistic
    S.bc_stat = T0-bias; % Bias-corrected test statistic
    S.SE = SE;           % Bootstrap standard error of the test statistic
    S.ci = ci.';         % Bootstrap confidence intervals of the test statistic
    S.prct = [m2,m1];    % Percentiles used to generate confidence intervals
    if any(diff(weights))
      S.weights = weights;
    else
      S.weights = [];
    end
    if isempty(clusters)
      S.strata = strata;
    else
      S.strata = [];
    end
    S.clusters = clusters;
    S.blocksize = blocksize;

  end

end

%--------------------------------------------------------------------------

function [T1, T2, U, idx] = boot1 (x, nboot, n, nvar, bootfun, T0, S, opt)

    % Extract required options structure fields
    weights = opt.weights;
    strata = opt.strata;
    blocksize = opt.blocksize;
    bandwidth = opt.bandwidth;
    R = opt.R;
    df = opt.df;
    stderr = opt.stderr;
    runmode = opt.runmode;
    
    % Initialize
    B = nboot(1);
    C = nboot(2);
    N = n*B;
    T1 = zeros(1,B);
    if C>0
      T2 = zeros(C,B);
      U = zeros(1,B);
    elseif ~isempty(stderr)
      T2 = [];
      U = zeros(1,B);
    else
      T2 = [];
      U = [];
    end
    X1 = cell(1,nvar);
    if nargout < 4
      idx = zeros(n,1);
    else
      idx = zeros(n,B);
    end

    % If applicable, prepare for stratified resampling
    if ~isempty(strata)
      % Get strata IDs
      gid = unique(strata);  % strata ID
      K = numel(gid);        % number of strata
      % Create strata matrix
      g = zeros(n,K);
      for k = 1:K
        g(:,k) = (strata == gid(k));
      end
      % Get strata sample and bootstrap sample set dimensions
      nk = sum(g).';   % strata sample sizes
      ck = cumsum(nk); % cumulative sum of strata sample sizes
      ik = [1;ck];     % strata boundaries
      Nk = nk*B;       % size of strata bootstrap sample set
      Ck = cumsum(Nk); % cumulative sum of strata bootstrap sample set sizes
    else
      ck = n;
      g = ones(n,1);
      K = 1;
      nk = n;
    end
    g = logical(g);

    % Prepare weights for resampling
    if any(diff(weights))
      if ~isempty(strata)
        % Calculate within-stratum weights
        c = zeros(n,1);
        for k = 1:K
          c = c + round(Nk(k) * g(:,k).*weights./sum(g(:,k).*weights));
          c(ik(k):ik(k+1),1) = cumsum(c(ik(k):ik(k+1),1));
          c(ik(k+1)) = Ck(k);
        end
      else
        % Calculate weights (no groups)
        c = cumsum(round(N * weights./sum(weights)));
        c(end) = N;
      end
      c = [c(1);diff(c)];
    else
      c = ones(n,1)*B;
    end

    % Perform bootstrap resampling
    try
      pool = gcp('nocreate');
    catch
      pool = [];
    end
    if isempty(pool)
      % Since first bootstrap is large, use a memory
      % efficient balanced resampling algorithm
      % If strata is provided, resampling is stratified
      for h = 1:B
        for i = 1:n
          k = sum(i>ck)+1;
          j = sum((rand(1) >= cumsum((g(:,k).*c)./sum(g(:,k).*c))))+1;
          if nargout < 4
            idx(i,1) = j;
          else
            idx(i,h) = j;
          end
          c(j) = c(j)-1;
        end
        for v = 1:nvar
          if nargout < 4
            X1{v} = x{v}(idx);
          else
            X1{v} = x{v}(idx(:,h));
          end
        end
        % Since second bootstrap is usually much smaller, perform rapid
        % balanced resampling by a permutation algorithm
        if C>0
          [U(h), T2(:,h)] = boot2 (X1, nboot, n, nvar, bootfun, T0, g, S, opt);
        end
        if ~isempty(bandwidth)
          % Apply smoothing using Gaussian/Student-t kernel
          noise = bsxfun(@times,mvtrnd(R,df,n),bandwidth);
          for v = 1:nvar
            X1{v} = X1{v} + noise(:,v);
          end
        end
        T1(h) = feval(bootfun,X1{:});
        if ~isempty(stderr)
          U(h) = stderr(X1{:});
        end
      end
    else
      % Matlab parallel mode
      % Perform ordinary resampling with replacement
      if nargout > 3
        error('No bootidx when operating ibootci in parallel mode')
      end
      % Prepare resampling weights 
      w = zeros(n,K);
      for k = 1:K
        w(:,k) = cumsum((g(:,k).*weights)./sum(g(:,k).*weights));
      end
      parfor h = 1:B
        idx = zeros(n,1);
        X1 = cell(1,nvar);
        for i = 1:n
          k = sum(i>ck)+1;
          j = sum((rand(1) >= w(:,k)))+1;
          idx(i,1) = j;
        end
        for v = 1:nvar
          X1{v} = x{v}(idx);
        end
        % Since second bootstrap is usually much smaller, perform rapid
        % balanced resampling by a permutation algorithm
        if C>0
          [U(h), T2(:,h)] = boot2 (X1, nboot, n, nvar, bootfun, T0, g, S, opt);
        end
        if ~isempty(bandwidth)
          % Apply smoothing using Gaussian/Student-t kernel
          noise = bsxfun(@times,mvtrnd(R,df,n),bandwidth);
          for v = 1:nvar
            X1{v} = X1{v} + noise(:,v);
          end
        end
        T1(h) = feval(bootfun,X1{:});
        if ~isempty(stderr)
          U(h) = stderr(X1{:});
        end
      end
    end

end

%--------------------------------------------------------------------------

function [U, T2] = boot2 (X1, nboot, n, nvar, bootfun, T0, g, S, opt)

    % Extract required options structure fields
    blocksize = opt.blocksize;
    runmode = opt.runmode;
    
    % Note that weights are not implemented here with iterated bootstrap

    % Prepare for block resampling (if applicable)
    if ~isempty(blocksize)
      x1 = cat_blocks(S.nvar,X1{:});
      blocksize = round(blocksize/2);
      X1 = split_blocks(x1,blocksize);
      nvar = S.nvar * blocksize;
      g = ones(n,1);
    end

    % Initialize
    C = nboot(2);

    % If applicable, prepare for stratified resampling
    K = size(g,2);    % number of strata
    nk = sum(g).';    % strata sample sizes
    ck = cumsum(nk);  % cumulative sum of strata sample sizes
    ik = [1;ck+1];    % strata boundaries (different definition to boot1)
    Nk = nk*C;        % size of strata bootstrap sample set

    % Rapid balanced resampling by permutation
    % If strata is provided, resampling is stratified
    idx = zeros(n,C);
    for k = 1:K
      tmp = (1:nk(k))'*ones(1,C);
      tmp = tmp(reshape(randperm(Nk(k)),nk(k),C));  % For compatibility with R2007
      tmp = tmp + ik(k) - 1;
      idx(ik(k): ik(k+1)-1,:) = tmp;
    end
    X2 = cell(1,nvar);
    for v = 1:nvar
      X2{v} = X1{v}(idx);
    end
    switch lower(runmode)
      case {'fast'}
        % Vectorized calculation of second bootstrap statistics
        T2 = feval(bootfun,X2{:});
      case {'slow'}
        % Calculation of second bootstrap statistics using a loop
        T2 = zeros(1,C);
        for i=1:C
          x2 = cellfun(@(X2)X2(:,i),X2,'UniformOutput',false);
          T2(i) = feval(bootfun,x2{:});
        end
    end
    U = interp_boot2(T2,T0,C);

end

%--------------------------------------------------------------------------

function  U = interp_boot2 (T2, T0, C)

    U = sum(T2<=T0);
    if U < 1
      U = 0;
    elseif U == C
      U = C;
    else
      % Quick linear interpolation to approximate asymptotic calibration
      t2 = zeros(1,2);
      I = (T2<=T0);
      if any(I)
        t2(1) = max(T2(I));
      else
        t2(1) = min(T2);
      end
      I = (T2>T0);
      if any(I)
        t2(2) = min(T2(I));
      else
        t2(2) = max(T2);
      end
      if (t2(2)-t2(1) == 0)
        U = t2(1);
      else
        U = ((t2(2)-T0)*U + (T0-t2(1))*(U+1)) /...
                (t2(2) - t2(1));
      end
    end
    U = U/C;

end

%--------------------------------------------------------------------------

function [m1, m2, S] = BC (B, func, x, T1, T0, alpha, S)

  % Note that alpha input argument is nominal coverage
    
  % Calculate bias correction z0
  z0 = norminv(sum(T1<T0)/B);

  % Calculate Bias-corrected percentiles
  z1 = norminv(0.5*(1+alpha));
  m1 = normcdf(2*z0+z1);
  z2 = norminv(0.5*(1-alpha));
  m2 = normcdf(2*z0+z2);
  S.z0 = z0;

end

%--------------------------------------------------------------------------

function [SSb, SSw, K, g, MSb, MSw, dk] = sse_calc (x, groups, nvar)

  % Calculate error components of groups

  % Initialize
  gid = unique(groups);  % group ID
  K = numel(gid);        % number of groups
  n = numel(x{1});
  g = zeros(n,K);
  bSQ = zeros(K,nvar);
  wSQ = zeros(n,nvar);
  center = zeros(K,nvar);
  % Calculate within and between group variances
  for k = 1:K
    % Create group matrix
    g(:,k) = (groups == gid(k));
    for v = 1:nvar
      center(k,v) = sum(g(:,k) .* x{v}) / sum(g(:,k));
      wSQ(:,v) = wSQ(:,v) + g(:,k).*(x{v}-center(k,v)).^2;
    end
  end
  for v = 1:nvar
    bSQ(:,v) = (center(:,v) - mean(center(:,v))).^2;
  end
  SSb = sum(bSQ);         % Between-group SSE
  SSw = sum(wSQ);         % Within-group SSE
  g = logical(g);         % Logical array defining groups

  % Calculate mean squared error (MSE) and representative cluster size
  if nargout > 4
    nk = sum(g).';
    MSb = (sum(nk.*bSQ))/(K-1);
    MSw = SSw/(n-K);
    dk = mean(nk) - sum((sum(g)-mean(nk)).^2)/((K-1)*sum(g(:)));
  end

end

%--------------------------------------------------------------------------

function [y, g] = unitmeans (x, clusters, nvar)

  % Calculate unit (cluster) means

  % Calculate number of levels of subsampling
  L = size(clusters,2);

  % Get IDs of unique clusters in lowest level
  gid = unique(clusters(:,L));
  K = numel(gid);

  % Initialize output variables
  g = zeros(K,L-1);
  y = cell(1,nvar);
  for v = 1:nvar
    y{v} = zeros(K,1);
  end

  % Calculate cluster means
  for k = 1:K

    % Find last level cluster members
    idx = find(clusters(:,L) == gid(k));

    % Compute cluster means
    for v = 1:nvar
      y{v}(k) = mean(x{v}(idx));
    end

    % Check data nesting
    if numel(unique(clusters(idx,L-1))) > 1
      error('Impossible hierarchical data structure')
    end

    % Redefine clusters
    g(k,:) = clusters(idx(1),1:L-1);

  end

end

%--------------------------------------------------------------------------

function [mu, Z, K, g] = clustmean (x, clusters, nvar)

  % Calculates shrunken cluster means and residuals for cluster bootstrap
  % See also bootclust function below

  % Center and scale data
  z = cell(1,nvar);
  for v = 1:nvar
    z{v} = (x{v} - mean(x{v})) / std(x{v});
  end

  % Calculate sum-of-squared error components
  [SSb, SSw, K, g] = sse_calc (z, clusters, nvar);
  SSb = sum(SSb);
  SSw = sum(SSw);

  % Calculate cluster means in the original sample
  mu = cell(1,nvar);
  for v = 1:nvar
    for k = 1:K
      mu{v}(k,:) = mean(x{v}(g(:,k),:));
    end
  end

  % Calculate shrunken cluster means from the original sample
  nk = sum(g).';
  dk = mean(nk) - sum((sum(g)-mean(nk)).^2)/((K-1)*sum(g(:)));
  c = 1 - sqrt(max(0,(K/(K-1)) - (SSw./(dk.*(dk-1).*SSb))));
  for v = 1:nvar
    mu{v} = bsxfun(@plus, c*mean(mu{v}),(1-c)*mu{v});
  end

  % Calculate residuals from the sample and cluster means
  Z = cell(1,nvar);
  for v = 1:nvar
    for k = 1:K
      Z{v}(g(:,k),:) = bsxfun(@minus, x{v}(g(:,k),:), mu{v}(k,:));
      Z{v}(g(:,k),:) = Z{v}(g(:,k),:) ./ sqrt(1-dk^-1);
    end
  end

end

%--------------------------------------------------------------------------

function T = bootclust (bootfun, K, g, runmode, mu, varargin)

  % Two-stage nonparametric bootstrap sampling with shrinkage
  % correction for clustered data [1-4].
  %
  % By resampling residuals, this bootstrap method can be used when
  % cluster sizes are unequal. However, cluster samples are assumed
  % to be taken from populations with equal variance. Not compatible
  % with bootstrap-t or bootstrap iteration.
  %
  % Reference5:
  %  [1] Davison and Hinkley (1997) Bootstrap Methods and their
  %       application. Chapter 3: pg 97-100
  %  [2] Ng, Grieve and Carpenter (2013) The Stata Journal.
  %       13(1): 141-164
  %  [3] Gomes et al. (2012) Medical Decision Making. 32(2): 350-361
  %  [4] Gomes et al. (2012) Health Econ. 21(9):1101-18
  
  % Calculate data dimensions
  Z = varargin{1};
  nvar = numel(Z);
  [n,reps] = size(Z{1});

  % Preallocate arrays
  bootmu = cell(1,nvar);
  X = cell(1,nvar);
  for v = 1:nvar
    X{v} = zeros(n,reps);
  end

  % Ordinary resampling with replacement of cluster means
  idx = ceil(K.*rand(K,reps));   % For compatibility with R2007
  for v = 1:nvar
    bootmu{v} = mu{v}(idx);
  end

  % Combine residuals with resampled cluster means
  for v = 1:nvar
    for k = 1:K
      X{v}(g(:,k),:) = bsxfun(@plus, Z{v}(g(:,k),:), bootmu{v}(k,:));
    end
  end

  % Calculate bootstrap statistic(s)
  switch lower(runmode)
    case {'fast'}
      T = feval(bootfun,X{:});
    case {'slow'}
      T = zeros(1,reps);
      for i = 1:reps
        x = cellfun(@(X) X(:,i),X,'UniformOutput',false);
        T(i) = feval(bootfun,x{:});
      end
  end

end

%--------------------------------------------------------------------------

function y = split_blocks (x, l)

  % Calculate data and block dimensions
  n = size(x{1},1);
  nvar = numel(x);

  % Create a matrix of circular, overlapping blocks
  % Ref: Politis and Romano (1991) Technical report No. 370
  y = cell(1,nvar);
  for v = 1:nvar
    y{v} = zeros(n,l);
    temp = cat(1,x{v},x{v}(1:l-1));
    for i = 1:n
      y{v}(i,:) = temp(i:i+l-1);
    end
  end
  y = cell2mat(y);
  y = num2cell(y,1);

end

%--------------------------------------------------------------------------

function y = cat_blocks (nvar, varargin)

  % Get data dimensions
  x = (varargin);
  N = numel(x);
  l = N/nvar;
  [n, reps] = size(x{1});

  % Concatenate blocks
  y = cell(1,nvar);
  for v = 1:nvar
    y{v} = zeros(N,reps);
    for i = 1:l
      y{v}(i:l:n*l,:) = x{(v-1)*l+i};
    end
    y{v} = y{v}(1:n,:);
  end

end

%--------------------------------------------------------------------------

function T = auxfun (bootfun, nvar, varargin)

  % Auxiliary function for block bootstrap
  X = varargin{1};
  Y = cat_blocks(nvar,X{:});
  T = bootfun(Y{:});

end

%--------------------------------------------------------------------------

function data = list2mat (varargin)

  % Convert comma-separated list input to matrix
  data = cell2mat(varargin);

end

%--------------------------------------------------------------------------

function [F, x] = empcdf (y, c)

  % Calculate empirical cumulative distribution function of y
  %
  % Set c to:
  %  1 to have a complete distribution with F ranging from 0 to 1
  %  0 to avoid duplicate values in x
  %
  % Unlike ecdf, empcdf uses a denominator of N+1

  % Check input argument
  if ~isa(y,'numeric')
    error('y must be numeric')
  end
  if all(size(y)>1)
    error('y must be a vector')
  end
  if size(y,2)>1
    y = y.';
  end

  % Create empirical CDF
  x = sort(y);
  N = sum(~isnan(y));
  [x,F] = unique(x,'rows','last');
  F = F/(N+1);

  % Apply option to complete the CDF
  if c > 0
    x = [x(1);x;x(end)];
    F = [0;F;1];
  end

end
