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
ttlem_params = default_ttlem_params();

%% lake params
lake_defs = [{'clipped_lost.tif', 3.523578530974025e+05, 1.612583134776515e+06}; ...
             {'clipped_jan.tif', 4.980370625e+05, 1.5489835e+06}];

core_depths = [1.38, 1.73]; % in m

[lakes, exp_volumes] = load_lakes_volumes(lake_defs, core_depths);
num_lakes = length(lakes);

%% run lakes
d_vals = linspace(0.02, 0.08, 2);
k_vals = linspace(0.000275, 0.001125, 2);

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
