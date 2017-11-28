function [params, output, marked_images] = thld_v5()
% THLD_V5 track horizontal line displacement and use shape language
% modeling to calculate shear  displacement, shear stress, shear strain,
% shear modulus, and shear energy dissipation in soft tissues
%
% Usage: 
%   [params, output, marked_image] = thld_v5();
%
% Notes:
% - Able to track multiple photobleached lines
% - Uses shape language modeling (least squares spline modeling using shape 
%   primitives) to smooth data & calculate derivatives (Relevant functions 
%   copyright (c) 2009, John D'Errico, available on the MATLAB file 
%   exchange. See below for link)
% - *** Before running this code, the user should modify the parameters in 
%   the subfunction "setup_params" according to the experiment. ***
%
% Author: 
% - Lena Bartell, March 2016. 
%   Based on previous versions by Mark Buckley (circa 2011-12?), Lena 
%   Bartell (circa 2013), and Corinne Henak (circa 2015)
%
% Relevant papers to cite when using this code:
% - John D'Errico. SLM - Shape Language Modeling. 15 Jun 2009 (Updated 29 
%   Apr 2014)
%   http://www.mathworks.com/matlabcentral/fileexchange/24443-slm-shape-language-modeling
% - Buckley et al. Localization of Viscous Behavior and Shear Energy 
%   Dissipation in Articular Cartilage Under Dynamic Shear Loading. J 
%   Biomech Eng. 2013 Mar; 135(3): 0310021–0310029.
% - Buckley et al. High-resolution spatial mapping of shear properties in 
%   cartilage. J Biomech. 2010. 43(4):796-800.

% Update History:
% - 2013-03-18, Rewrote Corinne's "THLD_multiple" code. Image import is now
%               faster and spline fits now use John D'Errico's Shape
%               Language Modeling. -Lena 
% - 2016-03-20, Corrected jumps in displacement phase that result from the
%               discontinuity of atan2 at +/-180 degrees. -Lena
% - 2106-04-05, Modify phase discontinuity fix. -Lena
% - 2016-04-11, Handled phase patching when there is no overlap. Modified
%               how phases are patched (more averaging over individual
%               lines first). Added constraint to force phase smoothing fit
%               to be zero at zero depth. -Lena 
% - 2016-04-20, Expanded knots parameters to allow for different regions
%               with different knot spacings. - Lena
% - 2016-04-28, Changed phase patching to remove outliers >4 STD away from
%               the mean. Added 'while' catch to ensure that the zero-depth
%               DC phase offset is not NAN. - Lena

% Display update to the screen
fprintf('THLD analysis... \n')
t1 = clock;

% Setup parameters 
% *** The user should modify these according to the experiment!
params = setup_params();

% Import the images, dividing them into tiles
[params, images] = import_tiles(params);

% Divide the tiles into strips and sum the data in each strip
[params, sum_data] = get_strips(params, images);

% Track the photobleached (PB) line(s) in each strip
[params, line_position] = track_lines(params, sum_data);

% USE THIS FOR DEBUGGING:
% Note: this is memory intensive, only use if you need to. To turn off,
% change "if 1" to "if 0"
if 0
    fprintf(' - Making movie of line position... ')
    t1 = clock;
    
    % Plot the line tracks over the image
    marked_images = plot_line_fit(params, images, line_position); 
    sz = size(marked_images);
    marked_images_squished = reshape(marked_images, [sz(1:3) sz(4)*sz(5)]);
    implay(marked_images_squished)
    
    t2 = clock;
    fprintf('Done (%2.2f s).\n', etime(t2,t1))
else
    marked_images = [];
end

% Fit the PB line position to extract its amplitude and phase
[params, amplitude, phase] = analyze_oscillation(params, line_position);

% Patch phases between tiles, so phase is continuous
[params, phase_patched] = patch_phases(params, phase);

% Ask user to remove any erroneous data
[params, full_disp_data] = remove_data(params, amplitude, phase_patched);

% Given all the raw displacement data, smooth it using shape-preserving 
% spline fits and calculate desired outputs (strain, shear modulus, etc)
[params, slmfit, data] = ...
    calculate_smooth_data(params, full_disp_data, amplitude, phase_patched);

% Combine processing results into a single output structure
output.sum_data = sum_data;
output.line_position = line_position;
output.amplitude = amplitude;
output.phase = phase;
output.phase_patched = phase_patched;
output.full_disp_data = full_disp_data;
output.slmfit = slmfit;
output.data = data;

% Plot some stuff before we go
plot_results(output)

% Say goodbye
t2 = clock;
fprintf('THLD analysis complete (%2.2f s total).\n', etime(t2,t1))

function params = setup_params()
% *** The user should modify these according to the experiment!

% Parameters for...
% - importing image
params.image.path = 'G:\Samples\Control 2\Control 2 (4-21-2016).tif'; % path to TIF stack that will be analyzed
params.image.umppx = 1.31; % Pixel size, in microns per pixel

% - dividing image into tiles
params.tile.number = 5; % total number of tiles (integer)
params.tile.frame_spacing = 1000; % number of frames between tiles in raw video (integer)
params.tile.analysis_start = 94; % for each tile, start analyzing at this frame (integer; best if the PB lines are straight in this frame & it is an integer number of cycles)
params.tile.analysis_length = 700; % for each tile, analyze this number of frames
params.tile.horizontal_jump.um = 587; % distance moved between each tile, Units: microns

% - masking background at articular surface & deepzone
params.mask.AS = 209; % number of columns to crop off the first tile (at AS), Units: px
params.mask.DZ = 390; % number of columns to crop off the last tile (at DZ), Units: px

% - finding photobleached (PB) lines
params.line.number = 5; % number of lines
params.line.first = 90;  % number of rows between the top of the image and the topmost PB line in first image, Units: px
params.line.spacing.um = 80; % distance between neighboring PB lines, Units: um
params.line.width.px = 8; % approximate width of each PB line, Units: px

% - analyzing oscillation
params.data.TDISconstant = 0.0080; % force constant for black TDIS (silver=0.0262), Last recalibrated Oct 20, 2013 (by CRH), Units: N/um
params.data.frequency =  1; % freqency of TDIS/tissue oscillation, Units: Hz
params.data.frame_rate = 41; % confocal data collection frame rate, Units: Frames per second
params.data.area = 16e-6; % Cross section area, Units: m^2
params.data.resolution.px = 10; % number of pixels to group into 1 depth data point (i.e. 1 strip), Units: px (integer)
params.data.spacing.px = 5; % number of pixels to move between each depth data point (i.e. each strip), Units: px (integer)

% - smoothing data & computing strain. 
% Note: smoothing is done via a "shape language model" (SLM), using the
% functions written by John D'Errico, available on the MATLAB file exchange
params.slm.degree = 3; % Degree of smoothing fit. Options: 0 (step function), 1, or 3
params.slm.monotonic = true; % If true, force amplitude spline fit to be monotonically decreasing.
params.slm.knots_step.um = [50 300]; %Units: um. Spacing between knots in spline fit, if this is a vector, divide the range into that many sections with the given spaicing in each
params.slm.knots_switch.um = [250]; %Units: um. Specify where to divide between different knot-spacing regions (this vector should be 1 element shorter than knots_step.um and should be evenly divisible by its associated step-size.
params.slm.show_plot = true; % If true, show diagnostic plots of the SLM fit (w/ data, knots, and resulting fit)

% - calculate other parameters based on the above info
params.tile.time_points = (1:params.tile.analysis_length)/...
    params.data.frame_rate; %time of each tile frame, Units: seconds
params.tile.horizontal_jump.px = params.tile.horizontal_jump.um / params.image.umppx;

params.data.resolution.px = round(params.data.resolution.px); % make sure resolution is an integer in px, Units: px
params.data.resolution.um = params.data.resolution.px * params.image.umppx; % Units: um
params.data.spacing.px = round(params.data.spacing.px); % make sure resolution is an integer in px, Units: px
params.data.spacing.um = params.data.spacing.px * params.image.umppx; % Units: um
params.data.frames_per_cycle = params.data.frame_rate / params.data.frequency; % Number of frames for one oscillation, Units: frames

params.line.spacing.px = params.line.spacing.um / params.image.umppx; % Units: px
params.line.width.um = params.line.width.px * params.image.umppx;
params.line.search.px = params.line.spacing.px * 0.5; % distance from estimate to look for the PB line, Units: px

function [params, images] = import_tiles(params)

fprintf(' - Reading images... ')
t1 = clock;

% gather image information
iminfo = imfinfo(params.image.path);
params.image.height = iminfo(1).Height;
params.image.width = iminfo(1).Width;
params.image.pages = length(iminfo);

% pull some variables from parameters
nTiles = params.tile.number;
tilespacing = params.tile.frame_spacing;
tilestart = params.tile.analysis_start;
tilelength = params.tile.analysis_length;

% establish which pages to read from the tif stack
tile_startix = ((1:nTiles)-1) * tilespacing + tilestart;
tile_endix = tile_startix + tilelength - 1;

% check that we aren't reading too many pages
if tile_endix(end)>params.image.pages
    error('Attempting to read too may pages from the input image')
end

% initialize output variable "images", which will hold the raw image data
% The dimensions of "images" correspond to: 
%   (image height, image width, analysis frame/page, tile)
firstim = imread(params.image.path, tile_startix(1), 'Info', iminfo);
images = zeros(params.image.height, params.image.width, ...
    tilelength, nTiles, 'like', firstim);
for tt = 1:nTiles
    
    pages = tile_startix(tt):tile_endix(tt);
    
    for pp = 1:tilelength
        images(:,:,pp,tt) = imread(params.image.path, pages(pp), ...
            'Info', iminfo);
    end
end

% flip horizontally so articular surface is on the left and deep zone is on
% the right
images = flipdim(images, 2);

t2 = clock;
fprintf('Done (%2.2f s).\n', etime(t2,t1))

function [params, sum_data] = get_strips(params, images)

fprintf(' - Gathering strips... ')
t1 = clock;

% - extract parameters
width = params.data.resolution.px;
spacing = params.data.spacing.px;

% - initialize outputs
[params.strip.start, params.strip.end, params.strip.center.px, ...
    params.strip.center.um, params.strip.depth.um, params.strip.depth.px, ...
    sum_data] = deal( cell(1, params.tile.number) );
params.strip.number = NaN(1, params.tile.number);

% Establish where each strip starts & ends in each tile & relative to all
% tiles
for tt = 1:params.tile.number
    
    % image column where strips start and end in each tile, Units: px
    firstcol = 1;
    lastcol = params.image.width;
    if tt==1
        firstcol = params.mask.AS;
    end
    if tt==params.tile.number
        lastcol = params.image.width - params.mask.DZ;
    end
    
    % column range for each strip in current tile
    params.strip.start{tt} = firstcol:spacing:(lastcol-width);
    params.strip.end{tt} = params.strip.start{tt} + width - 1;
    params.strip.center.px{tt} = mean([params.strip.start{tt}; params.strip.end{tt}]);
    params.strip.center.um{tt} = params.strip.center.px{tt} * params.image.umppx;
    params.strip.number(tt) = length(params.strip.start{tt});
    params.strip.depth.px{tt} = ...
        params.tile.horizontal_jump.px*(tt-1) + ...
        params.strip.center.px{tt} - ...
        params.mask.AS ;
    params.strip.depth.um{tt} = params.strip.depth.px{tt} * params.image.umppx;    
end

% For each strip in each tile, sum the data in the given range
for tt = 1:params.tile.number
    
    % dimensions of sum data: image height x strip x time 
    sum_data{tt} = NaN(params.image.height, params.strip.number(tt), ...
        params.tile.analysis_length);
    
    % sum each strip
    for ss = 1:params.strip.number(tt)
        a = params.strip.start{tt}(ss);
        b = params.strip.end{tt}(ss);
        sum_data{tt}(:,ss,:) = sum(images(:,a:b,:,tt), 2);
    end
end

t2 = clock;
fprintf('Done (%2.2f s).\n', etime(t2,t1))

function [params, line_position] = track_lines(params, sum_data)

fprintf(' - Tracking PB lines... ')
t1 = clock;

% parameters
nTiles = params.tile.number;
nFrames = params.tile.analysis_length;
nLines = params.line.number;
add_lines = (0:nLines-1) * round(params.line.spacing.px);

% initialize output
line_position = cell(1, nTiles);

% Loop over each tile
for tt = 1:nTiles
    
    nStrips = params.strip.number(tt);
    line_position{tt} = NaN(nLines, nStrips, nFrames);
    
    % Loop over each strip
    for ss = 1:nStrips
        
        % Loop over each frame
        for ff = 1:nFrames
            
            % Estimate where topmost PB line is
            estimate = NaN;
            wider_to_start = false;
            if ff==1 && ss==1
                % where user clicked on the PB line
                estimate = params.line.first;
                wider_to_start = true;

            elseif ff==1 && ss>1
                % where topmost line was in the prev strip
                estimate = line_position{tt}(1, ss-1, ff);

            elseif ff>1 && ss==1
                % where topmost line was in the prev frame
                estimate = line_position{tt}(1, ss, ff-1);

            elseif ff>1 && ss>1
                % average where it was in prev strip and prev frame
                estimate = mean([line_position{tt}(1, ss-1, ff),...
                    line_position{tt}(1, ss, ff-1) ]);
            end
            all_estimates = round(estimate) + add_lines;
            
            % Fit to extract where all PB lines are
            line_position{tt}(:,ss,ff) = find_minima(...
                params, sum_data{tt}(:,ss,ff), all_estimates, wider_to_start);
        end
    end
end

t2 = clock;
fprintf('Done (%2.2f s).\n', etime(t2,t1))

function yloc = find_minima(params, data, estimates, wider_to_start)

% initialize output variable
yloc = NaN(1, params.line.number);

% Loop over each line 
for ii = 1:params.line.number
    
    % Estimate of current line's position
    currest = estimates(ii);
    
    % Find local minima in the neighborhood of the estimate and then fit
    % a parabola around this minimum to extract the "true" minimum location
    if wider_to_start
        delta_min = round(params.line.spacing.px/2); 
    else
        delta_min = round(params.line.spacing.px/4); 
    end
    delta_fit = round(params.line.width.px); % window around data min for fitting
    linepos = find_fit_min(data, currest, delta_min, delta_fit);
    
    % Check if current position is close enough to the estimate. 
    % If not close enough, try a few different variations and check again:
    if abs(linepos-currest) > params.line.search.px 
        
        % - try a wider window for finding the minimum 
        delta_min = round(params.line.spacing.px/3);
        delta_fit = round(params.line.width.px);
        linepos = find_fit_min(data, currest, delta_min, delta_fit);        
        
        % Check again
        if abs(linepos-currest) > params.line.search.px,   
            
            % - try a wider window for fitting the minimum 
            delta_min = round(params.line.spacing.px/4);
            delta_fit = round(params.line.width.px*2);
            linepos = find_fit_min(data, currest, delta_min, delta_fit);        
            
            % Check again
            if abs(linepos-currest) > params.line.search.px,
                
                % - try wider windows for finding and fitting 
                delta_min = round(params.line.spacing.px/3);
                delta_fit = round(params.line.width.px*2);
                linepos = find_fit_min(data, currest, delta_min, delta_fit);        
                
                % Check again. 
                % If the above didn't work, throw a warning and set the 
                % line position to equal the estimate
                if abs(linepos-currest) > params.line.search.px, 
                    fprintf(['WARNING: Line %d was found outside the search range.' ...
                        'Setting it to be the estimate instead.\n'], ii)
                    plot(data, 'k.-')
                    hold on
                    plot(currest, data(round(currest)), 'r*')
                    plot(linepos, data(round(currest)), 'md')
                    title('there was a problem. execution paused.')
                    legend('data','estimate','line found here')
                    waitfor(gcf)
                    linepos = currest;
                end
            end
        end
    end
    
    % store the current line position
    yloc(ii) = linepos;
    
end

function out = find_fit_min(data, estimate, delta_min, delta_fit)

% limits x data to this range
xmin = 1;
xmax = length(data);
estimate = round(estimate);
delta_min = round(delta_min);
delta_fit = round(delta_fit);

% find minimum in the data near the current estimate
xdata = (estimate-delta_min):(estimate+delta_min);
xdata((xdata<xmin) | (xdata>xmax)) = []; % remove data that is out of range
ydata = data(xdata);
[~, ix] = min(ydata);
estimate = xdata(ix);

% fit neighborhood around minimum with a 2nd order polynomial and
% extract the "true" minimum
xdata = (estimate-delta_fit):(estimate+delta_fit);
xdata((xdata<xmin) | (xdata>xmax)) = []; % remove data that is out of range
ydata = data(xdata);
p = mypoly2fit(xdata, ydata); 
out = p(2) / (-2*p(1));
    
function p = mypoly2fit(x, y)
%based on polyfit, but without warnings or other orders
%assuming degree 2, and double-type inputs

x = x(:);
y = y(:);

if length(x)~=length(y)
    error(message('MATLAB:mypoly2fit:XYSizeMismatch'))
end

% Construct Vandermonde matrix.
V(:,3) = ones(length(x),1);
for j = 2:-1:1
   V(:,j) = x.*V(:,j+1);
end

% Solve least squares problem.
p = V\y;

function [params, amplitude, phase] = analyze_oscillation(params, line_position)

fprintf(' - Analyzing oscillation... ')
t1 = clock;

% parameters
nTiles = params.tile.number;
nFrames = params.tile.analysis_length;
nLines = params.line.number;
frequency = params.data.frequency;
FPS = params.data.frame_rate; % frames per second
time_pts = (1:nFrames)' / FPS;
fitX = [sin(2*pi*frequency*time_pts) ...
    cos(2*pi*frequency*time_pts) ...
    ones(size(time_pts))];

% initialize output
[amplitude.px, amplitude.um, phase] = deal(cell(1, nTiles));

% Loop over each tile
for tt = 1:nTiles
    
    nStrips = params.strip.number(tt);
    [amplitude.px{tt}, amplitude.um{tt}, phase{tt}] = deal(NaN(nLines, nStrips));
    
    % Loop over each strip
    for ss = 1:nStrips
        
        %Loop over each line
        for ii = 1:nLines
            
            % gather sum data
            data = line_position{tt}(ii, ss, :);
            data = permute(data, [3 1 2]);
            
            % fit data to extract phase and amplitude
            % NOTE: frequency is taken as known
            %       Fit model:
            %       y(t) = A*sin(2*pi*f*t + phi) + C
            %            = a*sin(2*pi*f*t) + b*cos(2*pi*f*t) + c
            K = fitX \ data;
            amplitude.px{tt}(ii, ss) = sqrt(K(1)^2 + K(2)^2);
            phase{tt}(ii, ss) = atan2(K(2), K(1)); %radians, in range: [-pi,pi]
            
            % USE THIS FOR DEBUGGING
            if 0
                % Plot oscillation fit result
                fity = fitX * K;
                figure('color', 'w')
                plot(time_pts, data, 'k.')
                hold on
                plot(time_pts, fity, '-')
                title(...
                    sprintf('Oscillation fit, tile %d, strip %d, line %d,',...
                    tt, ss, ii))
                waitfor(gcf)
            end
        end
    end
    
    % also convert from px to um
    amplitude.um{tt} = amplitude.px{tt} * params.image.umppx;
    
    % unwrap phase to correct for discontinuities near +/-180 deg
    phase{tt} = unwrap(phase{tt}, [], 2); % radians
end

t2 = clock;
fprintf('Done (%2.2f s).\n', etime(t2,t1))

function [params, phase] = patch_phases(params, phase)

fprintf(' - Patching phases... ')
t1 = clock;

% parameters
nTiles = params.tile.number;

% Loop over tiles
for tt = 1:nTiles
    
    % remove obvious outliers 
    data = phase{tt};
    tmp = ones(1, size(data, 2));
    avg_data = mean(data, 2) * tmp;
    std_data = std(data, [], 2) * tmp;
    keep = abs(data-avg_data) < 4*std_data;
    data(~keep) = NaN;
    phase{tt} = data;
    
end

% Loop over tiles 
for tt = 2:nTiles
    
    % phase data from each tile
    data1 = phase{tt-1};
    data2 = phase{tt};
    
    % find overlap between tiles
    d1 = params.strip.depth.px{tt-1};
    d2 = params.strip.depth.px{tt};
    nOverlap = sum(d1>=d2(1));
    
    % Check that there is enough overlap (5 pts minimum)
    nPtsMin = 5;
    if nOverlap < nPtsMin
        
        % if the overlap too small, use points at the edge of each region
        avg1 = nanmean(data1(:, end-nPtsMin+1:end), 2);
        avg2 = nanmean(data2(:, 1:nPtsMin), 2);
        
        % if the edges were all NaNs
        while ~any(isfinite(avg1)) || ~any(isfinite(avg2))
            % expand the region and try again
            nPtsMin = nPtsMin+1;            
            avg1 = nanmean(data1(:, end-nPtsMin+1:end), 2);
            avg2 = nanmean(data2(:, 1:nPtsMin), 2);
        end
    else
        
        % otherwise, find mean phase for full overlap region
        avg1 = nanmean(data1(:, d1>=d2(1)), 2);
        avg2 = nanmean(data2(:, d2<d1(end)), 2);
    end
    
    % correct for the difference
    tmp = ones(1, size(data2, 2));
    phase{tt} = phase{tt} - avg2*tmp + avg1*tmp;
end

% Subtract off the DC offset of each line
dc_offset = nanmean(cat(2, phase{:}), 2);
for tt = 1:nTiles
    tmp = ones(1, size(phase{tt}, 2));
    phase{tt} = phase{tt} - dc_offset*tmp;
end

% Loop over tiles and subtract the mean phase of the first few depth points
all_phases = cat(2, phase{:});
phase_d0 = nanmean(all_phases(:,1));
endix = 1;
while isnan(phase_d0)
    endix = endix+1;
    tmp = all_phases(:,1:endix);
    phase_d0 = nanmean(tmp(:));
end
for tt = 1:nTiles
    phase{tt} = phase{tt} - phase_d0;
end

t2 = clock;
fprintf('Done (%2.2f s).\n', etime(t2,t1))

function [params, data] = remove_data(params, A, phi)

fprintf(' - Removing erroneous data... ')
t1 = clock;

% gather all depth, amplitude & phase data
depth = cat(2, params.strip.depth.px{:});
depth = repmat(depth, [params.line.number 1]);
depth = depth(:);
amplitude = cat(2, A.px{:});
amplitude = amplitude(:);
phase = cat(2, phi{:});
phase = phase(:);
data1 = [depth amplitude];
data2 = [depth phase];

% brush to remove erroneous data
out1 = brush_data(data1);
out2 = brush_data(data2);

% collect output
data.amplitude.depth.px = out1(:,1);
data.amplitude.depth.um = data.amplitude.depth.px * params.image.umppx;
data.amplitude.px = out1(:,2);
data.amplitude.um = data.amplitude.px * params.image.umppx;
data.phase.depth.px = out2(:,1);
data.phase.depth.um = data.phase.depth.px * params.image.umppx;
data.phase.rad = out2(:,2);

t2 = clock;
fprintf('Done (%2.2f s).\n', etime(t2,t1))

function data = brush_data(data)
% Plot data, user brushes unwanted data, return wanted data only

% Plot all data vs depth, ask user to brush to remove any erroneous data
f = figure;
plot(data(:,1), data(:,2), 'k.')
ax = gca;
brush on
hmsg = msgbox('Highlight any erroneous data. Click OK when done');
waitfor(hmsg)

% Get brushed data
hchild = get(ax, 'Children');
brushed_data = NaN(0, 2);
for ii = 1:length(hchild)
    
    currh = hchild(ii).BrushHandles;
    
    if ~isempty(currh)
        
        currdata = currh.Children(1).VertexData;
        
        if ~isempty(currdata)
            brushed_data = cat(1, brushed_data, currdata(1:2,:)');
        end    
    end
end

% Remove brushed data
bix = ismember(data, brushed_data, 'rows');
data(bix,:) = [];

% Close figure
close(f)

function [params, slmfit, smoothdata] = calculate_smooth_data(params, data, tile_amplitude, tile_phi)

fprintf(' - Smoothing & calculating results... ')
t1 = clock;

% When reporting results, sample at these depths. Units: microns
mindepth = max( min(data.amplitude.depth.um), min(data.phase.depth.um) );
maxdepth = min( max(data.amplitude.depth.um), max(data.phase.depth.um) );
smoothdata.depth = mindepth : params.data.spacing.um : maxdepth; % units: um
smoothdata.UNITS.depth = 'microns';

% Gather data for the SLM spline fit for displacement amplitude
x = data.amplitude.depth.um;
y = data.amplitude.um;

% Create the knots for the SLM spline fit for displacement amplitude
% The knots vary based on the format of the input, potentially dividing
% the depth into different regions with different knot spacings.
nRegions = length(params.slm.knots_step.um);
if nRegions==1
    
    % specify the entire region at once
    knots = 0 : params.slm.knots_step.um : max(x);
    
    % force last knot to sit at the end of the data
    knots(end) = max(x); 
else
    
    % loop through each region and specify the knots
    knots = [];
    kstep = params.slm.knots_step.um;
    kstart = [0 params.slm.knots_switch.um(:)' max(x)];
    for rr = 1:nRegions
        knots = [knots, kstart(rr) : kstep(rr) : kstart(rr+1)-1];    
    end
    
    % force last knot to sit at the end of the data
    knots(end) = max(x);
end

% Create & sample the spline fit for displacement amplitude. Units: um
% Also calculate the first derivative d(u0)/d(z)
slmfit.disp_amplitude = slmengine(x, y, ...
    'degree', params.slm.degree, ...
    'knots', knots, ...
    'decreasing', ifout(params.slm.monotonic,'on','off'),...
    'endconditions', 'notaknot', ...
    'plot', ifout(params.slm.show_plot, 'on', 'off') );
smoothdata.displacement.amplitude = ...
    slmeval(smoothdata.depth, slmfit.disp_amplitude, 0);
smoothdata.UNITS.displacement.amplitude = 'microns';
smoothdata.displacement.amplitude_derivative = ...
    slmeval(smoothdata.depth, slmfit.disp_amplitude, 1);
smoothdata.UNITS.displacement.amplitude_derivative = 'none (micron per micron)';

% Gather data for the SLM spline fit for displacement phase
x = data.phase.depth.um;
y = data.phase.rad;

% Create the knots for the SLM spline fit for displacement phase
% The knots vary based on the format of the input, potentially dividing
% the depth into different regions with different knot spacings.
nRegions = length(params.slm.knots_step.um);
if nRegions==1
    
    % specify the entire region at once
    knots = 0 : params.slm.knots_step.um : max(x);
    
    % force last knot to sit at the end of the data
    knots(end) = max(x); 
else
    
    % loop through each region and specify the knots
    knots = [];
    kstep = params.slm.knots_step.um;
    kstart = [0 params.slm.knots_switch.um(:)' max(x)];
    for rr = 1:nRegions
        knots = [knots, kstart(rr) : kstep(rr) : kstart(rr+1)-1];    
    end
    
    % force last knot to sit at the end of the data
    knots(end) = max(x);
end

% Create & sample the spline fit for displacement phase. Units: rad
% Also calculate the first derivative d(phi)/d(z)
slmfit.disp_phase = slmengine(x, y, ...
    'degree', params.slm.degree, ...
    'knots', knots, ...
    'decreasing', 'off',...
    'leftvalue', 0,...
    'endconditions', 'notaknot', ...
    'plot', ifout(params.slm.show_plot, 'on', 'off') );
smoothdata.displacement.phase = ...
    slmeval(smoothdata.depth, slmfit.disp_phase, 0);
smoothdata.UNITS.displacement.phase = 'radians';
smoothdata.displacement.phase_derivative = ...
    slmeval(smoothdata.depth, slmfit.disp_phase, 1);
smoothdata.UNITS.displacement.phase_derivative = 'radians per micron';

% Use SLM fits & other input data to calculate...
% - Math used below is based on "Approach #3" from 
%   Buckley et al. J Biomech Eng. 2013 Mar; 135(3): 0310021–0310029.

% - strain amplitude & phase
u0 = smoothdata.displacement.amplitude;
du0_dz = smoothdata.displacement.amplitude_derivative;
phi = smoothdata.displacement.phase;
dphi_dz = smoothdata.displacement.phase_derivative;
grad_u0cosphi = - du0_dz .* cos(phi) + u0 .* sin(phi) .* dphi_dz;
grad_u0sinphi = - du0_dz .* sin(phi) - u0 .* cos(phi) .* dphi_dz;
smoothdata.strain.amplitude = sqrt(grad_u0cosphi.^2 + grad_u0sinphi.^2);
smoothdata.UNITS.strain.amplitude = 'none';
smoothdata.strain.phase = atan2(grad_u0sinphi, grad_u0cosphi);
smoothdata.UNITS.strain.phase = 'radians';

% - stress amplitude & phase
smoothdata.stress.amplitude = params.data.TDISconstant * ...
    mean(tile_amplitude.um{end}(:,end)) / params.data.area; 
smoothdata.stress.phase = mean(tile_phi{end}(:,end)); 
smoothdata.UNITS.stress.amplitude = 'Pa';
smoothdata.UNITS.stress.phase = 'radians';

% - shear modulus
smoothdata.shear_modulus = ...
    smoothdata.stress.amplitude ./ smoothdata.strain.amplitude;
smoothdata.UNITS.shear_modulus = 'Pa';

% - phase lag between stress and strain (i.e. viscosity)
smoothdata.phase_lag = smoothdata.stress.phase - smoothdata.strain.phase;
smoothdata.UNITS.phase_lag = 'radians';

% - energy dissipation per cycle per unit volume
smoothdata.energy_dissipated.raw = - pi * smoothdata.stress.amplitude^2 * ...
    sin(smoothdata.phase_lag) ./ smoothdata.shear_modulus; % multiply by (-1) because we are looking at dissipation
smoothdata.UNITS.energy_dissipated.raw = 'Joules per meter^3 per cycle';
Ediss_total = sum(smoothdata.energy_dissipated.raw);
smoothdata.energy_dissipated.normalized = ...
    smoothdata.energy_dissipated.raw / Ediss_total;
smoothdata.UNITS.energy_dissipated.normalized = 'local fraction';

t2 = clock;
fprintf('Done (%2.2f s).\n', etime(t2,t1))

function out = ifout(cond,a,b)
% IFOUT Ternary operator
%
% OUT = IFOUT(COND,A,B) implements a ternary operator.
% if cond == true then out = a
% if cond == false then out = b

if cond
    out = a ;
else
    out = b ;
end

function [] = plot_results(output)
% Plot results 
figure('color', 'w')

subplot(221)
x = output.full_disp_data.amplitude.depth.um;
y = output.full_disp_data.amplitude.um;
plot(x, y, 'k.')
x = output.data.depth;
y = output.data.displacement.amplitude;
hold on
plot(x, y, '-', 'linewidth', 2)
xlim([0 x(end)])
xlabel(sprintf('depth (%s)', output.data.UNITS.depth))
ylabel(sprintf('displacement amplitude (%s)', output.data.UNITS.displacement.amplitude))
legend({'raw data', 'SLM spline fit'})

subplot(222)
x = output.full_disp_data.phase.depth.um;
y = output.full_disp_data.phase.rad;
plot(x, y, 'k.')
x = output.data.depth;
y = output.data.displacement.phase;
hold on
plot(x, y, '-', 'linewidth', 2)
xlim([0 x(end)])
xlabel(sprintf('depth (%s)', output.data.UNITS.depth))
ylabel(sprintf('displacement phase (%s)', output.data.UNITS.displacement.phase))
legend({'raw data', 'SLM spline fit'})

subplot(223)
x = output.data.depth;
y = output.data.shear_modulus;
semilogy(x, y, '-', 'linewidth', 2)
xlim([0 x(end)])
xlabel(sprintf('depth (%s)', output.data.UNITS.depth))
ylabel(sprintf('shear modulus (%s)', output.data.UNITS.shear_modulus))

subplot(224)
x = output.data.depth;
y = output.data.energy_dissipated.normalized;
stairs(x, y, 'linewidth', 2)
xlim([0 x(end)])
ylim([0 1])
xlabel(sprintf('depth (%s)', output.data.UNITS.depth))
ylabel(sprintf(['energy dissipated per volume' char(10) 'per cycle (%s)'], output.data.UNITS.energy_dissipated.normalized))

function marked_images = plot_line_fit(params, images, line_position)
% show where the line was found in each image

marked_images = repmat(permute(images, [1 2 5 3 4]), [1 1 3]);

nLines = params.line.number;
nFrames = params.tile.analysis_length;
nTiles = params.tile.number;
imsize = [params.image.height params.image.width 3 nFrames nTiles];
colors = lines(nLines)*255;

for tt = 1:nTiles
    
    % get subscript in image of line positions
    x = repmat(params.strip.center.px{tt}', [1 nFrames]);
    x = round(x);
    nStrips = size(x,1);
    [~, frame] = ndgrid(1, 1:nFrames);
    frame = repmat(frame, [nStrips 1]);
    c = ones(size(x));
    tile = tt*c;
    
    for ll = 1:nLines
        
        % continue to get subscript in image of line positions
        y = permute(line_position{tt}(ll,:,:), [2 3 1]);
        y = round(y);
        
        % convert to index, for each color        
        R = c*1;
        G = c*2;
        B = c*3;        
        Rix = sub2ind(imsize, y, x, R, frame, tile);
        Gix = sub2ind(imsize, y, x, G, frame, tile);
        Bix = sub2ind(imsize, y, x, B, frame, tile);
        
        % set color for each pixel of the line
        marked_images(Rix) = colors(ll, 1);
        marked_images(Bix) = colors(ll, 2);
        marked_images(Gix) = colors(ll, 3);
    end
end











