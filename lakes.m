

%{
given
- DEM landscape file
- basin outlet

calculate filters:
- basin
- lake


1. apply gaussian filter to DEM file. why?
2. fill sinks


. calculate drainage basin
. extrapolate lake

}

%}


%{
  results: basin filter, lake filter, expected volume
%}

%
% function [sed] = calculate_sediment(dem, outlet_index, )
%
%
%
%
% end

close all, clear all, clc


project_dir = '/Users/pboone/Dropbox/Projects/Classes/GEOL 2049/project-lake-sediment-analysis/';

lake.dem = GRIDobj([project_dir, 'clipped_lost.tif']);
lake.dem_filled = NaN;
lake.outlet.x = 3.523578530974025e+05;
lake.outlet.y = 1.612583134776515e+06;
lake.outlet.z = NaN; % calculate after gaussian filter is applied, in case it changes
lake.outlet.index = coord2ind(lake.dem, lake.outlet.x, lake.outlet.y);
lake.sediment_depth = 1.38;

%% prep dem file GRIDobj

% apply gaussian filter
hsize = 10;
sigma = 2;
h = fspecial('gaussian', hsize, sigma);
lake.dem.Z=imfilter(lake.dem.Z, h, 'replicate');
lake.outlet.z = lake.dem.Z(lake.outlet.index);

% fill sinks
lake.dem_filled = fillsinks(lake.dem);

topo.cell_area = lake.dem.cellsize^2;

topo.flow.direction = FLOWobj(lake.dem_filled);
% topo.flow.acceleration = flowacc(topo.flow.direction)
topo.gradient = gradient8(lake.dem);

%% determine lake size

lake_gradient_threshold = 0.03;
lake_depth = 1;
d = drainagebasins(topo.flow.direction, lake.outlet.index);
topo.drainage_basin_filter = d.Z;
topo.lake_filter = lake.dem.Z*0;
topo.lake_filter(topo.drainage_basin_filter == 1 & ...
                 topo.gradient.Z < lake_gradient_threshold & ...
                 lake.dem.Z < lake.outlet.z + lake_depth ...
                 )=1;

% calculate sediment

num_lake_cells = length(find(topo.lake_filter == 1));
num_basin_cells = length(find(topo.drainage_basin_filter == 1)) - num_lake_cells;

%% calculate sediment volume from lake size

sediment_volume = lake.sediment_depth * num_lake_cells * topo.cell_area;
disp(sediment_volume);
