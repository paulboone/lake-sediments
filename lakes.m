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

addpath ~/workspace/export_fig
%}

close all, clear all, clc

%% TTLEM vars
ttlem_params.TimeSpan = 6000;
ttlem_params.TimeStep = 6000;

% Hillslope processes, parameters values
ttlem_params.m = 0.5;
ttlem_params.n = 1;

% Fixed or variable drainage network through time
ttlem_params.DrainDir = 'variable';

ttlem_params.AreaThresh = 2e5;     % channel contributing area threshold

% Threshold slopes
ttlem_params.Sc = 1;
ttlem_params.Sc_unit = 'tangent';

% Output
ttlem_params.ploteach=inf;
ttlem_params.saveeach=1;

%% run lake
project_dir = '/Users/pboone/Dropbox/Projects/Classes/GEOL 2049/project-lake-sediment-analysis/';
lost_lake = Lake();
lost_lake.load_from_geotiff([project_dir, 'clipped_lost.tif'], 3.523578530974025e+05, 1.612583134776515e+06);

jan_lake = Lake();
jan_lake.load_from_geotiff([project_dir, 'clipped_jan1.tif'], 4.980370625e+05, 1.5489835e+06);

v_exp_jan = jan_lake.calculate_sediment_volume_from_core(1.73);
v_exp_lost = lost_lake.calculate_sediment_volume_from_core(1.38);

d_vals = linspace(0, 0.3, 10);
k_vals = linspace(0, 0.01, 10);

for i=1:length(d_vals);
  disp(['Progress: ' num2str(i-1) '/' num2str(length(d_vals))]);
  for j=1:length(k_vals)
    ttlem_params.D = d_vals(i);
    ttlem_params.Kw = k_vals(j);

    ttlem_params = ttlemset(ttlem_params);
    v_mod = jan_lake.calculate_sediment_volume_via_model(ttlem_params);
    chi2(i,j) = ((v_mod - v_exp_jan)^2)/(v_exp_jan^2);

    v_mod = lost_lake.calculate_sediment_volume_via_model(ttlem_params);
    chi2(i,j) = chi2(i,j) + ((v_mod - v_exp_lost)^2)/(v_exp_lost^2);

    close all;
  end
end

imagesc(k_vals, d_vals, log10(chi2))
colorbar
shg

% still need to verify units on these:
% xlabel('K [m^2yr^{-1}]');
% ylabel('D [m^2yr^{-1}]');
