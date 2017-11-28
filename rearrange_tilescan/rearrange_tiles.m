function tilescan = rearrange_tiles(raw_tilescan, num_tiles_vertical, num_tiles_horizontal)
% tilescan = REARRANGE_TILES(raw_tilescan, num_tiles_vertical, num_tiles_horizontal)
% Re-arrange the tiles in RAW_TILESCAN so they are in the correct order.
%
% Inputs:
%   raw_tilescan            Image data (NxMx1 or NxMx3 matrix) with
%                           mis-arranged tiles
%
%   num_tiles_vertical      Integer number of tiles in the vertical (y)
%                           direction in RAW_TILESCAN
%
%   num_tiles_horizontal    Integer number of tiles in the horizontal (x)
%                           direction in RAW_TILESCAN
% 
% Outputs:
%   tilescan                Fixed (i.e. re-arranged) tilescan
%
% Notes:
% - If the tiles are arranged as:
%        1  2  3  4
%        5  6  7  8
%        9 10 11 12
%   then this function re-arranges them to be arranged as:
%       12  9  6  3 
%       11  8  5  2
%       10  7  4  1
%
% Example implementation:
% tilescan = imread('tilescan.tif');
% nY = 2;
% nX = 2;
% fixed_tilescan = rearrange_tiles(tilescan, nY, nX);
%
% Lena Bartell, Nov 2017

% divide raw image into tiles
nY = num_tiles_vertical;
nX = num_tiles_horizontal;
tiles = mat2cell(raw_tilescan, ...
    repmat(size(raw_tilescan, 1)/nY, 1, nY),...
    repmat(size(raw_tilescan, 2)/nX, 1, nX),...
    repmat(size(raw_tilescan, 3)/1, 1, 1));

% fix order
tiles = flip(flip(permute(tiles, [2 1]), 1), 2);

% convert back to full image
tilescan = cell2mat(tiles);

