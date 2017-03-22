function [sediment_volume, final_dem] = calculate_sediment_volume_via_model(lake, ttlem_params)
  ttlem_params = ttlemset(ttlem_params);

  % no uplift lem
  uplift_dem = lake.dem;
  uplift_dem.Z = uplift_dem.Z * 0;
  ttlem_results = ttlem(lake.dem, uplift_dem, ttlem_params);
  final_dem = ttlem_results.H1;

  % find landscape Z differences
  dzdt = lake.dem.Z - final_dem.Z;
  dzdt(lake.drainage_basin_filter ~= 1 | lake.lake_filter == 1) = NaN;
  sediment_volume = nansum(dzdt(:)) * lake.cell_area;
end