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

%% test for auto locate outlet sill
% lake_defs = [{'clipped_lost.tif', NaN, NaN}; ...
%              {'clipped_jan.tif', NaN, NaN}];


lake_defs = [{'clipped_lost.tif', 3.523578530974025e+05, 1.612583134776515e+06}; ...
             {'clipped_jan.tif', 4.980370625e+05, 1.5489835e+06}];

core_depths = [1.38, 1.73]; % in m

d_vals = linspace(0.07, 0.08, 9);
k_vals = linspace(0.000875, 0.001125, 9);

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
chi2all = NaN(length(d_vals), length(k_vals));
for i=1:length(d_vals)
  disp(['Progress: ' num2str(i-1) '/' num2str(length(d_vals))]);
  for j=1:length(k_vals)
    [chi2(:,i,j), chi2all(i,j)] = calc_chi_for_lakes(lakes, exp_volumes, ttlem_params, k_vals(j), d_vals(i));
  end
end

%% plot
for i=1:num_lakes
  lake = lakes(i);
  lakechi2 = squeeze(chi2(i,:,:));
  plotsave_chi2(lake.lake_name, lakechi2, k_vals, d_vals);
end

plotsave_chi2('all', chi2all, k_vals, d_vals);
