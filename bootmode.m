% Function file: [H, P, h] = bootmode (x, m, B, kernel)
%
% This function tests whether the distribution underlying the univariate
% data in vector x has m modes. The method employs the smooth bootstrap
% as described [1].
%
% The parsimonious approach is to consider a successively increasing
% number of modes until the null hypothesis (H0) is accepted (i.e. H=0),
% where H0 corresponds to the number of modes being equal to m.
%
% x is the vector of data
%
% m is the number of modes for hypothesis testing
%
% B is the number of bootstrap replicates
%
% kernel can be 'Gaussian' (default) or 'Epanechnikov'
%
% H=0 indicates that the null hypothesis cannot be rejected at the 5%
% significance level.  H=1 indicates that the null hypothesis can be
% rejected at the 5% level.
%
% P is the achieved significance level using the bootstrap test.
%
% h is the critical bandwidth (i.e. the smallest bandwidth achievable to
% obtain a kernel density estimate with m modes)
%
% Bootstrap iteration is not implemented for this test.
%
%    Bibliography:
%    [1] Efron and Tibshirani. Chapter 16 Hypothesis testing with the
%         bootstrap in An introduction to the bootstrap (CRC Press, 1994)
%
%  bootmode v1.1.1.0 (22/04/2019)
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


%% Stamp data example used in reference [1]
%% From stamp dataset in bootstrap R package
%x=[0.060;0.064;0.064;0.065;0.066;0.068;0.069;0.069;0.069;0.069;0.069;0.069;0.069;0.070;0.070;0.070;
%0.070;0.070;0.070;0.070;0.070;0.070;0.070;0.070;0.070;0.070;0.070;0.070;0.070;0.070;0.070;0.070;
%0.070;0.070;0.070;0.070;0.070;0.070;0.070;0.071;0.071;0.071;0.071;0.071;0.071;0.071;0.071;0.071;
%0.071;0.071;0.071;0.071;0.071;0.071;0.071;0.071;0.071;0.071;0.071;0.072;0.072;0.072;0.072;0.072;
%0.072;0.072;0.072;0.072;0.072;0.072;0.072;0.072;0.072;0.072;0.072;0.072;0.072;0.072;0.072;0.072;
%0.072;0.072;0.072;0.072;0.072;0.072;0.072;0.072;0.072;0.072;0.072;0.073;0.073;0.073;0.073;0.073;
%0.073;0.073;0.073;0.073;0.073;0.073;0.074;0.074;0.074;0.074;0.074;0.074;0.074;0.074;0.074;0.074;
%0.075;0.075;0.075;0.075;0.075;0.075;0.075;0.075;0.075;0.075;0.075;0.075;0.075;0.075;0.075;0.075;
%0.075;0.075;0.075;0.075;0.076;0.076;0.076;0.076;0.076;0.076;0.076;0.076;0.076;0.076;0.076;0.076;
%0.076;0.076;0.076;0.076;0.076;0.076;0.077;0.077;0.077;0.077;0.077;0.077;0.077;0.077;0.077;0.077;
%0.077;0.078;0.078;0.078;0.078;0.078;0.078;0.078;0.078;0.078;0.078;0.078;0.078;0.078;0.078;0.078;
%0.078;0.078;0.078;0.078;0.078;0.078;0.078;0.078;0.079;0.079;0.079;0.079;0.079;0.079;0.079;0.079;
%0.079;0.079;0.079;0.079;0.079;0.079;0.079;0.079;0.079;0.079;0.079;0.079;0.079;0.079;0.079;0.079;
%0.079;0.079;0.079;0.079;0.079;0.079;0.079;0.079;0.079;0.079;0.079;0.079;0.079;0.079;0.079;0.079;
%0.079;0.079;0.080;0.080;0.080;0.080;0.080;0.080;0.080;0.080;0.080;0.080;0.080;0.080;0.080;0.080;
%0.080;0.080;0.080;0.080;0.080;0.080;0.080;0.080;0.080;0.080;0.080;0.080;0.080;0.080;0.080;0.080;
%0.080;0.080;0.080;0.080;0.080;0.080;0.080;0.081;0.081;0.081;0.081;0.081;0.081;0.081;0.081;0.081;
%0.081;0.081;0.081;0.081;0.081;0.081;0.082;0.082;0.082;0.082;0.082;0.082;0.082;0.082;0.082;0.082;
%0.082;0.082;0.082;0.082;0.082;0.082;0.082;0.082;0.083;0.083;0.083;0.083;0.083;0.083;0.083;0.084;
%0.084;0.084;0.085;0.085;0.086;0.086;0.087;0.088;0.088;0.089;0.089;0.089;0.089;0.089;0.089;0.089;
%0.089;0.089;0.089;0.090;0.090;0.090;0.090;0.090;0.090;0.090;0.090;0.090;0.091;0.091;0.091;0.092;
%0.092;0.092;0.092;0.092;0.093;0.093;0.093;0.093;0.093;0.093;0.094;0.094;0.094;0.095;0.095;0.096;
%0.096;0.096;0.097;0.097;0.097;0.097;0.097;0.097;0.097;0.098;0.098;0.098;0.098;0.098;0.099;0.099;
%0.099;0.099;0.099;0.100;0.100;0.100;0.100;0.100;0.100;0.100;0.100;0.100;0.100;0.100;0.100;0.100;
%0.100;0.100;0.101;0.101;0.101;0.101;0.101;0.101;0.101;0.101;0.101;0.102;0.102;0.102;0.102;0.102;
%0.102;0.102;0.102;0.103;0.103;0.103;0.103;0.103;0.103;0.103;0.104;0.104;0.105;0.105;0.105;0.105;
%0.105;0.106;0.106;0.106;0.106;0.107;0.107;0.107;0.108;0.108;0.108;0.108;0.108;0.108;0.108;0.109;
%0.109;0.109;0.109;0.109;0.109;0.109;0.110;0.110;0.110;0.110;0.110;0.110;0.110;0.110;0.110;0.110;
%0.110;0.111;0.111;0.111;0.111;0.112;0.112;0.112;0.112;0.112;0.114;0.114;0.114;0.115;0.115;0.115;
%0.117;0.119;0.119;0.119;0.119;0.120;0.120;0.120;0.121;0.122;0.122;0.123;0.123;0.125;0.125;0.128;
%0.129;0.129;0.129;0.130;0.131]

% Results from multimodality testing using this
% function with the Gaussian kernel compare quantitatively with reference [1]
% Inference is similar when using an Epanechnikov kernel
%
% [H, P, h] = bootmode(x,m,2000,'Gaussian')
% m	 h       P     H
% 1  0.0068  0.00  1
% 2  0.0032  0.32  0 <- smallest m where H0 accepted
% 3  0.0030  0.07  0
% 4  0.0028  0.00  1
% 5  0.0026  0.00  1
% 6  0.0024  0.00  1
% 7  0.0015  0.46  0
%
% [H, P, h] = bootmode(x,m,2000,'Epanechnikov')
% m	 h       P     H
% 1  0.0169  0.02  1
% 2  0.0069  0.57  0 <- smallest m where H0 accepted
% 3  0.0061  0.52  0
% 4  0.0057  0.41  0
% 5  0.0056  0.19  0
% 6  0.0052  0.17  0
% 7  0.0052  0.06  0

function [H, P, h] = bootmode (x, m, B, kernel)

  if nargin < 2
    error('A minimum of 2 input arguments are required')
  end
  if nargin < 3
    B = '2000';
  end
  if nargin < 4
    kernel = 'Gaussian';
  end

  % Vectorize the data
  x = x(:);
  n = numel(x);

  % Find critical bandwidth
  [criticalBandwidth] = findCriticalBandwidth (x, m, kernel);
  h = criticalBandwidth;

  % Random resampling with replacement from a smooth estimate of the distribution
  idx = randi(n,n,B);
  Y = x(idx);
  xvar = var(x,1); % calculate sample variance
  Ymean = ones(n,1) * mean(Y);
  X = Ymean + (1 + h^2/xvar)^(-0.5) * (Y - Ymean + (h * randn(n,B)));

  % Calculate bootstrap statistic of the bootstrap samples
  % Boostrap statistic is the number of modes (faster)
  f = zeros(200,B);
  for j = 1:B
    [f(:,j),t] = kde(X(:,j),h, kernel);
  end
  % Compute number of modes in the kernel density estimate of the bootstrap samples
  mboot = sum(diff(sign(diff(f)))<0);
  % Approximate achieved significance level (ASL) from bootstrapping number of modes
  P = sum(mboot > m)/B;

  % Accept or reject null hypothesis
  if P > 0.05
    H = 0;
  elseif P < 0.05
    H = 1;
  end

end

%--------------------------------------------------------------------------

function [criticalBandwidth] = findCriticalBandwidth (x, mref, kernel)

  if mref > numel(x)
    error('the number of modes cannot exceed the number of datapoints')
  end

  % Calculate starting value for bandwidth
  % The algorithm sets the initial bandwidth so that half
  % of the sorted, unique data points are well separated
  xsort = sort(x);
  xdiff = diff(xsort);
  h = median(xdiff(xdiff>0))/(2*sqrt(2*log(2)));

  m = inf;
  while m > mref
    % Increase h
    h = h * 1.01; % Increase query bandwidth by 1%
    [y, t] = kde (x, h, kernel);
    m = sum(diff(sign(diff(y)))<0);
  end
  criticalBandwidth = h;

end

%--------------------------------------------------------------------------

function [f, t] = kde (x, h, kernel)

  % Vectorize the data
  x = x(:);

  % Default properties of t
  n = numel(x);
  tmin = min(x) - 3 * h; % define lower limit of kernel density estimate
  tmax = max(x) + 3 * h; % define upper limit of kernel density estimate

  % Kernel density estimator
  t = linspace(tmin,tmax,200)';
  f = zeros(200,1);
  for i = 1:n
    xi = ones(200,1) * x(i);
    u = (t-xi)/h;
    if strcmpi(kernel,'Epanechnikov')
      % Epanechnikov (parabolic) kernel
      K = @(u) max(0,(3/4)*(1-u.^2));
    elseif strcmpi(kernel,'Gaussian')
      % Gaussian kernel
      K = @(u) 1/sqrt(2*pi)*exp(-0.5*u.^2);
    end
    f = f + K(u);
  end
  f = f * 1/(n*h);

end
