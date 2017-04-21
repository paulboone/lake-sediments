

function lake = locate_outlet(geotiff_path, lake_point_x, lake_point_y)
  lake = Lake();
  lake.load_from_geotiff(geotiff_path, lake_point_x, lake_point_y);

  figure('Color',[1,1,1]);
  imagesc(lake.dem_filled);

  title(['Lake footprint, drainage basin and sill for ', lake.lake_name]);
  xlabel(['Lake cells: ', num2str(lake.num_lake_cells), char(10), ...
          'Lake basin cells: ', num2str(lake.num_basin_cells - lake.num_lake_cells), char(10), ...
          'Basin / lake Ratio: ', num2str((lake.num_basin_cells - lake.num_lake_cells)/lake.num_lake_cells)]);

  hold on;

  green1 = cat(3, zeros(size(lake.dem.Z)), 0.75*ones(size(lake.dem.Z)), zeros(size(lake.dem.Z)));
  h1 = imshow(green1, 'XData', lake.dem.georef.SpatialRef.XWorldLimits, ...
                     'YData', flip(lake.dem.georef.SpatialRef.YWorldLimits));


  green2 = cat(3, zeros(size(lake.dem.Z)), 0.5*ones(size(lake.dem.Z)), zeros(size(lake.dem.Z)));
  h2 = imshow(green2, 'XData', lake.dem.georef.SpatialRef.XWorldLimits, ...
                     'YData', flip(lake.dem.georef.SpatialRef.YWorldLimits));
  plot(lake.outlet.x,lake.outlet.y, 'r+','LineWidth',5);

  hold off;

  set(gca, 'YDir', 'normal');
  set(gca, 'Visible', 'on');

  set(h1, 'AlphaData', single(lake.drainage_basin_filter))
  set(h2, 'AlphaData', lake.lake_filter)
end
