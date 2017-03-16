

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
lake.outlet.index = coord2ind(lake.dem, lake.outlet.x, lake.outlet.y)
lake.sediment_depth = 1.38;

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

sediment_volume = lake.sediment_depth * num_lake_cells * topo.cell_area;


% Al=length(find(lk==1));
% % Aw=A.Z(IX)-Al
% Aw=nansum(B.Z(B.Z==1 & lk~=1))
% pa=DEM.cellsize*DEM.cellsize;
% %%
% D=Dv(i);
% Vemp(i)=D*Al*pa;
% Vemp_u(i)=Vemp(i)/2; %the case of tringular cs rather than rect



%
% %%flow direction object
% FD = FLOWobj(DEMf);
%
% %%flow acc
% A = flowacc(FD);
% G=gradient8(DEM);
%
% %%compute drainage basin for this outlet
% x=xv(i);
% y=yv(i);
% IX=coord2ind(DEM,x,y);%get index of outlet
%
% zo=DEM.Z(IX);%outlet elev
% B=drainagebasins(FD,IX);%get a binary grid of basin cells;
% Bv{i}=B;
%
% %%detect lake
% lk=DEM.Z*0;
% lk(B.Z==1 & G.Z<0.03 & DEM.Z<zo+1)=1;
% lkv{i}=lk;
% %
% Zb=DEM.Z*NaN;
% Zb(B.Z==1)=DEM.Z(B.Z==1);
% %

% imagesc(Zb);
% axis equal
% shg
% %     pause
%
% %%vol of sed in time frame
% Al=length(find(lk==1));
% % Aw=A.Z(IX)-Al
% Aw=nansum(B.Z(B.Z==1 & lk~=1))
% pa=DEM.cellsize*DEM.cellsize;
% %%
% D=Dv(i);
% Vemp(i)=D*Al*pa;
% Vemp_u(i)=Vemp(i)/2; %the case of tringular cs rather than rect












% imagesc(lake.dem)
% colorbar
% disp(lake.dem.Z(lake.outlet.index))
% dem1 = lake.dem
% % initial_colormap = colormap
%
% pause
%
% figure
%
%
% % DEMs.Z=imfilter(DEM.Z, h, 'replicate');
% % DEM=DEMs;
% disp(lake.dem.Z(lake.outlet.index))
% imagesc(lake.dem)
% % colormap(initial_colormap);
% colorbar
%
% dem2 = lake.dem


% DEM=DEMs;
