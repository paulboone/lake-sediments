

function lake = locate_outlet(geotiff_path)
  lake = Lake();
  lake.load_from_geotiff(geotiff_path, NaN, NaN);

  figure;
  imagesc(lake.dem_filled);

  hold on;

  green = cat(3, zeros(size(lake.dem.Z)), 0.5*ones(size(lake.dem.Z)), zeros(size(lake.dem.Z)));
  h = imshow(green, 'XData', lake.dem.georef.SpatialRef.XWorldLimits, ...
                    'YData', flip(lake.dem.georef.SpatialRef.YWorldLimits));
  plot(lake.outlet.x,lake.outlet.y, 'r+','LineWidth',5);

  hold off;

  set(gca, 'YDir', 'normal');
  set(gca, 'Visible', 'on');

  set(h, 'AlphaData', lake.lake_filter) 
end
