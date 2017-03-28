close all, clear all, clc;

project_dir = '/Users/pboone/Dropbox/Projects/Classes/GEOL 2049/project-lake-sediment-analysis/examples/';
l = Lake();
l.load_from_geotiff([project_dir, 'clipped_jan.tif'], NaN, NaN);

figure;
imagesc(l.dem_filled);

hold on;

green = cat(3, zeros(size(l.dem.Z)), 0.5*ones(size(l.dem.Z)), zeros(size(l.dem.Z)));
h = imshow(green, 'XData', l.dem.georef.SpatialRef.XWorldLimits, 'YData', flip(l.dem.georef.SpatialRef.YWorldLimits));
plot(l.outlet.x,l.outlet.y, 'r+','LineWidth',5);

hold off;

set(gca, 'YDir', 'normal');
set(gca, 'Visible', 'on');

set(h, 'AlphaData', l.lake_filter) 
