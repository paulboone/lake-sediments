

function lake = init_lake(geotiff_path, outlet_x, outlet_y)

  % constants for refactoring into a class
  gaussian_hsize = 10;
  gaussian_sigma = 2;
  lake_gradient_threshold = 0.03;
  lake_depth = 1;

  %% lake vars
  lake.dem = GRIDobj(geotiff_path);
  lake.dem_filled = NaN;
  lake.cell_area = lake.dem.cellsize^2;
  lake.outlet.x = outlet_x;
  lake.outlet.y = outlet_y;
  lake.outlet.z = NaN; % calculate after gaussian filter is applied, in case it changes
  lake.outlet.index = coord2ind(lake.dem, lake.outlet.x, lake.outlet.y);

  lake.flow_direction = NaN;
  lake.flow_acceleration = NaN;
  lake.gradient = NaN;
  lake.drainage_basin_filter = NaN;
  lake.lake_filter = NaN;

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

  lake.num_lake_cells = length(find(lake.lake_filter == 1));
  %num_basin_cells = length(find(topo.drainage_basin_filter == 1)) - num_lake_cells;
end
