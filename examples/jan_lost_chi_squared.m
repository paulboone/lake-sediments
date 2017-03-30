%{
Requires topotoolbox to be on the path.

Last used with version 2.1-108-g0f8ab25, from commit 0f8ab2575bae1fb88f927ad1e3f6f141987bfa58

Should be something like:
addpath ~/workspace/topotoolbox/
addpath ~/workspace/topotoolbox/utilities
addpath ~/workspace/topotoolbox/topoapp
addpath ~/workspace/topotoolbox/DEMdata
addpath ~/workspace/topotoolbox/ttlem
addpath ~/workspace/topotoolbox/ttlem/LateralDisplacement

addpath '~/Dropbox/Projects/Classes/GEOL 2049/project-lake-sediment-analysis/'


addpath ~/workspace/export_fig
%}

close all, clear all, clc

%% TTLEM vars
ttlem_params.TimeSpan = 6000;
ttlem_params.TimeStep = 6000;
ttlem_params.m = 0.5;
ttlem_params.n = 1;
ttlem_params.DrainDir = 'variable';
ttlem_params.AreaThresh = 2e5;     % channel contributing area threshold
ttlem_params.Sc = 1;
ttlem_params.Sc_unit = 'tangent';
ttlem_params.ploteach=inf;
ttlem_params.saveeach=1;

%% lake params
lake_defs = [{'clipped_lost.tif', 3.523578530974025e+05, 1.612583134776515e+06}; ...
             {'clipped_jan.tif', 4.980370625e+05, 1.5489835e+06}];

core_depths = [1.38, 1.73]; % in m

d_vals = linspace(0.0, 0.3, 10);
k_vals = linspace(0.0, 0.01, 10);

%% run lakes
num_lakes = size(lake_defs,1);

lakes(1,num_lakes) = Lake();
exp_volumes = zeros(num_lakes);
for i = 1:num_lakes
  lake_def = lake_defs(i,:);
  lakes(i).load_from_geotiff(lake_def{:});
  exp_volumes(i) = lakes(i).calculate_sediment_volume_from_core(core_depths(i));
end

chi2 = NaN(num_lakes, length(d_vals), length(k_vals));

for i=1:length(d_vals)
  disp(['Progress: ' num2str(i-1) '/' num2str(length(d_vals))]);
  for j=1:length(k_vals)
    ttlem_params.D = d_vals(i);
    ttlem_params.Kw = k_vals(j);

    ttlem_params = ttlemset(ttlem_params);

    for l=1:num_lakes
      v_mod = lakes(l).calculate_sediment_volume_via_model(ttlem_params);
      chi2(l,i,j) = ((v_mod - exp_volumes(l))^2)/(exp_volumes(l)^2);
    end

    close all;
  end
end


%% plot
for i=1:num_lakes
  lake = lakes(i);
  lakechi2 = squeeze(chi2(i,:,:));
  plotsave_chi2(lake.lake_name, lakechi2, k_vals, d_vals);
end

function plotsave_chi2(lake_name, lakechi2, k_vals, d_vals)
  figure
  imagesc(k_vals, d_vals, log10(lakechi2))
  title(['Chi^2 graph of varying K,D values (', lake_name, ')'])
  xlabel('K [m^2yr^{-1}]');
  ylabel('D [m^2yr^{-1}]');
  xticks(k_vals);
  yticks(d_vals);
  caxis([-6, 0])
  colorbar
  shg
  file_basename = [lake_name,'_k', num2str(k_vals(1)), '_', num2str(k_vals(end)), '_d', num2str(d_vals(1)), '_', num2str(d_vals(end))];
  save([file_basename, '.mat'], 'k_vals', 'd_vals', 'lakechi2');
  export_fig(['chi2_', file_basename,  '.png']);
end
