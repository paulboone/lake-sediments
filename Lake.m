

classdef Lake < handle
  properties

    gaussian_hsize
    gaussian_sigma
    lake_gradient_threshold
    lake_depth_threshold
    lake_name

  end

  properties (SetAccess = private, GetAccess = public)
    dem
    dem_filled
    cell_area
    outlet
    flow_direction
    flow_acceleration
    gradient
    drainage_basin_filter

    lake_point
    lake_max_depth
    lake_filter


    num_lake_cells
  end

  methods
    function obj = Lake()
      % set defaults:
      obj.gaussian_hsize = 10;
      obj.gaussian_sigma = 2;
      obj.lake_gradient_threshold = 0.03;
      obj.lake_depth_threshold = 1;
    end

    function load_from_geotiff(obj, geotiff_path, outlet_x, outlet_y)
      if isempty(obj.lake_name)
        [~,obj.lake_name,~] = fileparts(geotiff_path);
        obj.lake_name = strrep(obj.lake_name, '_', '-');
      end
      %% lake vars
      obj.dem = GRIDobj(geotiff_path);
      obj.cell_area = obj.dem.cellsize^2;

      % apply gaussian filter
      h = fspecial('gaussian', obj.gaussian_hsize, obj.gaussian_sigma);
      obj.dem.Z=imfilter(obj.dem.Z, h, 'replicate');

      % fill sinks
      obj.dem_filled = fillsinks(obj.dem);
      obj.flow_direction = FLOWobj(obj.dem_filled);
      obj.gradient = gradient8(obj.dem);

      % find point in lake and max depth: lake should include point with deepest fill
      z_diff = obj.dem.Z - obj.dem_filled.Z;
      [~, obj.lake_point.index] = min(z_diff(:));
      [obj.lake_point.x,obj.lake_point.y] = ind2coord(obj.dem, obj.lake_point.index);
      obj.lake_max_depth = z_diff(obj.lake_point.index);

      % set or determine outlet
      if ~isnan(outlet_x) && ~isnan(outlet_y)
        obj.outlet.x = outlet_x;
        obj.outlet.y = outlet_y;
        obj.outlet.index = coord2ind(obj.dem, obj.outlet.x, obj.outlet.y);
      else
        % automatically determine outlet by finding deepest sink and
        % following the flow pathway. outlet is first point in stream
        % that has lower elevation
        flowpath = obj.flow_direction.flowpathextract(obj.lake_point.index);
        outlet_stream = flowpath(obj.dem_filled.Z(flowpath) < obj.dem_filled.Z(flowpath(1)));
        obj.outlet.index = outlet_stream(1);
        [obj.outlet.x, obj.outlet.y] = ind2coord(obj.dem, obj.outlet.index);
      end

      obj.outlet.z = obj.dem.Z(obj.outlet.index);

      %% determine lake size
      d = drainagebasins(obj.flow_direction, obj.outlet.index);
      obj.drainage_basin_filter = d.Z;
      obj.lake_filter = obj.dem.Z*0;
      obj.lake_filter(obj.drainage_basin_filter == 1 & ...
                       obj.gradient.Z < obj.lake_gradient_threshold & ...
                       obj.dem.Z < obj.outlet.z + obj.lake_depth_threshold ...
                       )=1;
      %% gradient condition commented out for now: REVISIT
      %
      obj.num_lake_cells = length(find(obj.lake_filter == 1));
      % num_basin_cells = length(find(topo.drainage_basin_filter == 1)) - num_lake_cells;
    end


    function sediment_volume = calculate_sediment_volume_from_core(obj, sediment_depth)
      sediment_volume = sediment_depth * obj.num_lake_cells * obj.cell_area;
    end

    function [sediment_volume, final_dem] = calculate_sediment_volume_via_model(obj, ttlem_params)
      ttlem_params = ttlemset(ttlem_params);

      % no uplift lem
      uplift_dem = obj.dem;
      uplift_dem.Z = uplift_dem.Z * 0;
      ttlem_results = ttlem(obj.dem, uplift_dem, ttlem_params);
      final_dem = ttlem_results.H1;

      % find landscape Z differences
      dzdt = obj.dem.Z - final_dem.Z;
      dzdt(obj.drainage_basin_filter ~= 1 | obj.lake_filter == 1) = NaN;
      sediment_volume = nansum(dzdt(:)) * obj.cell_area;
    end
  end
end
