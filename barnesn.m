function [Vq, params, Xq] = barnesn(X, V, Xv, varargin)
% BARNESN Barnes smoothing interpolation of unstructured data
%   Vq = BARNESN(X, V, Xv) returns the smoothing interpolation of
%   D-dimensional observations V(X) at query points Xq. Query points Xq are
%   created by meshing the vectors in the cell array Xv that define the
%   grid in each dimension. Smoothing interpolation is performed using
%   the Koch form of Barnes objective analysis [2]. Roughly, (in 2D) the
%   interpolated value (vq) at gridpoint (xq, yq) is determined as a
%   weighted-sum of the values (v) at data points (x, y), based on the
%   gaussian weighting function exp(-r^2 / s / g^j), where r is the
%   euclidian distance from (xq, yq) to (x, y), s is the Gaussian Variance,
%   and g is the Convergence Parameter. 
%
%   Required inputs:
%   X       M x D array of data observation points, where M is the number
%           of observations and D is the dimensionality. 
%   Y       1 x D vector specifying the data value observed at each
%           data observation point in X  
%   Xv      D-element cell array of grid vectors {xgv1, xgv2,...xgvD} where
%           xgvi specifies the sample grid points along dimension i.
%
%   Vq = BARNESN(..., Name,Value) additionally specifies optional
%   analysis parameters using Name/Value pairs. In most cases, the default
%   parameter values are best. Optional Name/Value pairs include:  
% 
%   'Iterations'            Positive integer >= 2. The number of iterations
%                           to perform. The default is 2. This default
%                           value enables the fastest optimal convergence.
%                           Thus, this default value should be used in most
%                           cases.
% 
%   'ConvergenceParameter'  Positive double in the range (0 1]. The default
%                           value, 0.2. Practically it should be in the
%                           range (0.2, 1), where smaller ensures faster
%                           convergence. Thus, this default value should be
%                           used in most cases.
%
%   'GaussianVariance'      Positive double. The variance of the smoothing
%                           gaussian weights, s. The default value is
%                           optimized based on both the data spacing and
%                           the ConvergenceParameter. Thus, this default
%                           value should be used in most cases.
% 
%   'Quiet'                 Boolean. If true, suppress interation output.
%                           Default is false.
%
%   [Vq, params, Xq] = BARNESN(...) also returns the parameters
%   used for analysis in the structure params, as well as the full list of
%   query grid points in each dimension, stored in the cell array Xq, where
%   Xq{i} contains a D-dimensional array of grid coordinates along the i-th
%   dimension. Note that Xq grid points conform to the format produced by
%   "NDGRID" not "MESHGRID".
%
%   See also griddata ndgrid scatteredInterpolant interp2 interp3
%
%   Bibliography:     
%   [1] Barnes, Stanley L. "Mesoscale objective map analysis using weighted
%       time-series observations." (1973)
%   [2] Koch, Steven E., Mary DesJardins, and Paul J. Kocin. "An
%       interactive Barnes objective map analysis scheme for use with
%       satellite and conventional data." Journal of Climate and Applied
%       Meteorology 22.9 (1983): 1487-1503.
%   [3] Daley, Roger. Atmospheric data anlysis. No. 2. Cambridge University
%       Press, 1993.
%
%   Required toolboxes:
%       Statistics and Machine Learning Toolbox
%
%   Tested on: 
%       MATLAB 2016a (version 9.0)
%       with Statistics and Machine Learning Toolbox version 10.2
%
%   Inspired by Pierce, Stephen. "Barnes Objective Analysis."  MATLAB File
%   Exhange, File ID #28666. 

%% Example implementation:
%
%     % Data observation points
%     num_pts = 500;
%     num_dim = 2;
%     X = rand(num_pts, num_dim);
% 
%     % Sampling grid vector and associated query points
%     Xv = repmat({0:0.02:1}, [1, num_dim]);
%     Xq0 = cell(num_dim, 1);
%     [Xq0{:}] = ndgrid(Xv{:});
%     Xq = cellfun(@(a){a(:)}, Xq0);
%     Xq = cat(2, Xq{:});
% 
%     % Underlying function is f(x_i) = x_1 + x_2 + ...
%     % Based on observation and grid points, calculate:
%     % - True value at grid points
%     % - True value plus gaussian noise at observation points (i.e. generated data)
%     % - Smoothing interpolation at grid points, based on generated data
%     Vq0 = reshape(sum(Xq, 2), size(Xq0{1}));
%     V = sum(X,2) + randn(num_pts,1)/40;
%     [Vq, params, Xq1] = barnesn(X, V, Xv);
% 
%     % Plot the results
%     switch num_dim
%         case 1
%             fig = figure('color', 'w');
%             ax = gca;
%             ax.Box = 'on';
%             ax.NextPlot = 'add';
%             hpts = plot(ax, X, V, 'k.', 'DisplayName', 'Observation points (with noise)');
%             p0 = plot(ax, Xq, Vq0(:), 'c-', 'DisplayName', 'Underlying function');
%             p = plot(ax, Xq, Vq(:), 'm.-', 'DisplayName', 'Smoothed interpolation');
%             lg = legend([hpts, p0, p]);
% 
%         case 2
%             % Plot results as a surface
%             fig = figure('color', 'w');
%             ax = gca;
%             ax.View = [20 20];
%             ax.Box = 'on';
%             ax.NextPlot = 'add';
%             hpts = plot3(ax, X(:,1), X(:,2), V, 'k.', 'DisplayName', 'Observation points (with noise)');
%             p0 = surf(Xq0{:}, Vq0, 'DisplayName', 'Underlying function');
%             p0.FaceColor = 'cyan';
%             p0.FaceAlpha = 0.5;
%             p = surf(Xq0{:}, Vq, 'DisplayName', 'Smoothed interpolation');
%             p.FaceColor  = 'magenta';
%             p.FaceAlpha = 0.5;
%             lg = legend([hpts, p0, p]);
% 
%         case 3
%             % Plot results as an isosurface
%             fig = figure('color', 'w');
%             ax = gca;
%             ax.View = [20 20];
%             ax.Box = 'on';
%             ax.NextPlot = 'add';
%             iso_pt = 1.5;
%             hpts = plot3(ax, X(:,1), X(:,2), X(:,3), 'k.', 'DisplayName', 'Observation points (with noise)');
%             p0 = patch(isosurface(Xq0{[2 1 3]}, Vq0, iso_pt), 'DisplayName', ...
%                 sprintf('Underlying function isosurface at %0.1f', iso_pt));
%             p0.FaceColor = 'cyan';
%             p0.FaceAlpha = 0.5;
%             p = patch(isosurface(Xq0{[2 1 3]}, Vq, iso_pt), 'DisplayName', ...
%                 sprintf('Smoothed interpolation isosurface at %0.1f', iso_pt));
%             p.FaceColor  = 'magenta';
%             p.FaceAlpha = 0.5;
%             lg = legend([hpts, p0, p]);
%     end

%% Setup for analysis
% Parse inputs
[X, V, Xq, params] = parse_inputs(X, V, Xv, varargin{:});

% Setup gaussian weight factors for both iterations and normalize them
r = pdist2(Xq, X); % distance from each grid pt to each data pt
outer_data = ones(1, params.nData);
W = cell(2, 1);
for ii = 0:params.Iterations-1
    w = exp( -r.^2 / params.GaussianVariance / params.ConvergenceParameter^ii );
    sum_w = sum(w, 2) * outer_data;
    W{ii+1} = w ./ sum_w;
end

%% First pass
ii = 0;
outer_grid = ones(params.nGrid, 1);
f = outer_grid * V';
Vq = sum(W{ii+1}.*f, 2);

%% Subsequent passes

% Calculate error at observation points
Xcell = mat2cell(X, params.nData, ones(params.D, 1));
[Verr, params] = calc_Verror(V, Xv, Vq, Xcell, params, ii);

% For each additional iteration...
for ii = 1:params.Iterations-1

    % apply analysis to error field and 
    f = outer_grid * Verr';

    % Add this error field to the result from the previous iteration
    Vq = Vq + sum(W{ii+1}.*f, 2);

    % Interpolate and calculate new error at observation points
    [Verr, params] = calc_Verror(V, Xv, Vq, Xcell, params, ii);
    
end

%% Reshape outputs before returning
Vq = reshape(Vq, params.grid_size);
Xq = mat2cell(Xq, [params.nGrid], ones(1,params.D));
for dd = 1:params.D
    Xq{dd} = reshape(Xq{dd}, params.grid_size);
end

function [X, V, Xq, params] = parse_inputs(varargin)
% Parse & prep inputs and parameters

% Parse required inputs 
p = inputParser;
addRequired(p, 'X') % data points (observation station locations)
addRequired(p, 'V') % data values (station observations)
addRequired(p, 'Xv') % grid vectors
addParameter(p, 'GaussianVariance', NaN) % variance of gaussian weights
addParameter(p, 'ConvergenceParameter', 0.2, @(x) x>0 && x<=1) % convergence parameter, Default is 0.2 (based on Koch 1983)
addParameter(p, 'Iterations', 2, @(x) x>=2) % Number of iterations. Default is 2 (based on Koch 1983)
addParameter(p, 'Quiet', false); % print reports to the screen
parse(p, varargin{:})
X = p.Results.X;
V = p.Results.V(:);
Xv = p.Results.Xv(:);
params.GaussianVariance = p.Results.GaussianVariance;
params.ConvergenceParameter = p.Results.ConvergenceParameter;
params.Iterations = p.Results.Iterations;
params.Quiet = p.Results.Quiet;

% Check that X and V have the same number of data points
if size(X,1) ~= size(V,1);
    error('MATLAB:InputParser:ArgumentFailedValidation', ...
        ['The sizes of ''V'' and ''X'' do not match. ' ...
        'V must have one value for each row in X']);
end

% Remove data points with any Nans/Infs
keep = all(isfinite(X), 2);
X = X(keep,:);
V = V(keep);

% Setup parameters, store variable sizes
params.D = size(X,2);  % number of dimensions 
params.nData = size(X,1); % number of data observation points
params.grid_size = cellfun(@numel, Xv(:)'); % size of grid in each dimension
if params.D==1
    params.grid_size(2) = 1;
end
params.nGrid = prod(params.grid_size); % number of total query points
params.RMSE_avg = NaN(params.Iterations, 1); % store average RMS error of each iteration

% Setup the query grid points as a Q-by-D array, where Q is the number of
% query grid points and D is the number of dimensions
Xq = cell(params.D,1);
[Xq{:}] = ndgrid(Xv{:}); 
Xq = cellfun(@(a){a(:)}, Xq);
Xq = [Xq{:}];

% Estimate the average spacing between data points. If the data sampling is
% very non-uniform, this measure is better than (and much larger than) the
% average nearest-neighbor spacing. 
A = prod(max(X) - min(X)); % data extent "area"
M = params.nData;
params.data_spacing = sqrt(A)*(1+sqrt(M))/(M-1);

% Use data spacing to set optimal gaussian standard deviation, check if
% gaussian std is too small
optimal_var = (2*params.data_spacing/pi)^2 * params.ConvergenceParameter^(-params.Iterations);

if isnan(params.GaussianVariance)
    params.GaussianVariance = optimal_var;
end
if params.GaussianVariance < optimal_var
    warning('Gaussian variance is small. The optimal value is %f', optimal_var)
end
params.GaussianStd = sqrt(params.GaussianVariance);

% Warn if the grid spacing is not appropriate for the data spacing
% Limits from Koch 1983: 1/3 <= (dn/[grid spacing]) <= 1/2, where dn is the
% average nearest-neighbor spacing of the data points
min_grid_spacing = min(cellfun(@(a)min(diff(a)), Xv));
max_grid_spacing = max(cellfun(@(a)max(diff(a)), Xv));
if (min_grid_spacing/params.data_spacing) < (1/3) && ~params.Quiet
    warning('Grid spacing should be larger than [data spacing]/3 = %f. Smallest grid spacing: %f, data spacing: %f', ...
        params.data_spacing/3, min_grid_spacing, params.data_spacing)
end
if (max_grid_spacing/params.data_spacing) > (1/2) && ~params.Quiet
    warning('Note that grid spacing can be smaller than [data spacing]/2 = %f. Largest grid spacing: %f, data spacing: %f', ...
        params.data_spacing/2, max_grid_spacing, params.data_spacing)
end

% Warn if some grid points are far from the data points
r = pdist2(Xq, X); % distance from each grid pt to each data pt
if any(sum(r<=2*params.GaussianStd,2) < 3) && ~params.Quiet
    warning('Some grid points are far from any data points. Consider modifying the grid.');
end

function [Verr, params] = calc_Verror(V, Xv, Vq, Xcell, params, ii)
% Forward interpolate analyzed grid values to the observation locations and
% calculate the error between observed and interpolated values (at
% iteration ii)

% Interpolate and calculate error
Verr = V - interpn(Xv{:}, reshape(Vq, params.grid_size), Xcell{:});
params.RMSE_avg(ii+1) = nanmean(Verr.^2);

% Replace nans (from data points that are outside the grid) with zeros
outside = isnan(Verr);
Verr(outside) = 0;

% Print the average error to the screen
if ~params.Quiet
    fprintf('Barnes iteration %d/%d, average RMS error %0.5f\n', ...
        ii, params.Iterations-1, params.RMSE_avg(ii+1));
end


















