

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

% add topotoolbox paths - double check what is needed!
addpath ~/workspace/topotoolbox/
addpath ~/workspace/topotoolbox/utilities
addpath ~/workspace/topotoolbox/topoapp
addpath ~/workspace/topotoolbox/DEMdata
addpath ~/workspace/topotoolbox/ttlem
addpath ~/workspace/topotoolbox/ttlem/LateralDisplacement
addpath ~/workspace/export_fig
addpath ~/workspace/freezeColors


close all, clear all, clc




%{
load(geotiff_path, outlet_x, outlet_y, param_struct)

calculate_sediment_volume(sediment_depth)

%}


project_dir = '/Users/pboone/Dropbox/Projects/Classes/GEOL 2049/project-lake-sediment-analysis/';

gaussian_hsize = 10;
gaussian_sigma = 2;

lake_gradient_threshold = 0.03;
lake_depth = 1;

%% lake vars
lake.dem = GRIDobj([project_dir, 'clipped_lost.tif']);
lake.dem_filled = NaN;
lake.cell_area = lake.dem.cellsize^2;
lake.outlet.x = 3.523578530974025e+05;
lake.outlet.y = 1.612583134776515e+06;
lake.outlet.z = NaN; % calculate after gaussian filter is applied, in case it changes
lake.outlet.index = coord2ind(lake.dem, lake.outlet.x, lake.outlet.y);
lake.sediment_depth = 1.38;
lake.flow_direction = NaN;
% lake.flow_acceleration = NaN;
lake.gradient = NaN;
lake.drainage_basin_filter = NaN;
lake.lake_filter = NaN;


%% TTLEM vars
% Temporal domain
ttlem_params.TimeSpan = 6000;
ttlem_params.TimeStep = 6000;

% Hillslope processes, parameters values
ttlem_params.m = 0.5;
ttlem_params.n = 1;

% Fixed or variable drainage network through time
ttlem_params.DrainDir = 'variable';

% Threshold slopes
ttlem_params.Sc = 1;
ttlem_params.Sc_unit = 'tangent';

% Output
ttlem_params.ploteach=inf;
ttlem_params.saveeach=1;

ttlem_params.AreaThresh = 2e5;     % channel contributing area threshold, m�


%% PREP: INIT
%% prep dem file GRIDobj

% apply gaussian filter
h = fspecial('gaussian', gaussian_hsize, gaussian_sigma);
lake.dem.Z=imfilter(lake.dem.Z, h, 'replicate');
lake.outlet.z = lake.dem.Z(lake.outlet.index);

% fill sinks
lake.dem_filled = fillsinks(lake.dem);
lake.flow_direction = FLOWobj(lake.dem_filled);
lake.gradient = gradient8(lake.dem);
% lake.flow_acceleration = flowacc(lake.flow_direction)

%% determine lake size
d = drainagebasins(lake.flow_direction, lake.outlet.index);
lake.drainage_basin_filter = d.Z;
lake.lake_filter = lake.dem.Z*0;
lake.lake_filter(lake.drainage_basin_filter == 1 & ...
                 lake.gradient.Z < lake_gradient_threshold & ...
                 lake.dem.Z < lake.outlet.z + lake_depth ...
                 )=1;

num_lake_cells = length(find(lake.lake_filter == 1));
%num_basin_cells = length(find(topo.drainage_basin_filter == 1)) - num_lake_cells;

%% sediment_volume(sediment_depth)
%% calculate sediment volume from lake size

sediment_volume = lake.sediment_depth * num_lake_cells * lake.cell_area;
disp(sediment_volume);


%% TTLEM(lake, ttlem_params)

ttlem_params.D = 0.1;
ttlem_params.Kw = 0.05;

ttlem_params = ttlemset(ttlem_params);

% no uplift
uplift_dem = lake.dem;
uplift_dem.Z = uplift_dem.Z * 0;

ttlem_results = ttlem(lake.dem, uplift_dem, ttlem_params);
final_dem = ttlem_results.H1;

% find landscape Z differences
dzdt = lake.dem.Z - final_dem.Z;
dzdt(lake.drainage_basin_filter ~= 1 & lake.lake_filter == 1) = NaN;
volume = nansum(dzdt(:)) * lake.cell_area;
disp(volume);













%%
