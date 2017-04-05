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

global lakes ttlem_params exp_volumes

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

%% load lakes
num_lakes = size(lake_defs,1);

lakes = Lake();
lakes(1,num_lakes) = Lake();
exp_volumes = zeros(num_lakes);
for i = 1:num_lakes
  lake_def = lake_defs(i,:);
  lakes(i).load_from_geotiff(lake_def{:});
  exp_volumes(i) = lakes(i).calculate_sediment_volume_from_core(core_depths(i));
end

xmin = fminsearch(@sediment_volume_kd, [0.0012, 0.067 ]); 

disp(xmin);

function chi2 = sediment_volume_kd(kd)
  global lakes ttlem_params exp_volumes
  
  k = kd(1);
  d = kd(2);
  
  if k < 0 || d < 0
    chi2 = 1000;
  else
    [~, chi2] = calc_chi_for_lakes(lakes, exp_volumes, ttlem_params, k, d);
  end
  
  disp([k, d, chi2]);
    
  close all;
end

