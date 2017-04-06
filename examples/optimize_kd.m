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
ttlem_params = default_ttlem_params();

%% lake params
lake_defs = [{'clipped_lost.tif', 3.523578530974025e+05, 1.612583134776515e+06}; ...
             {'clipped_jan.tif', 4.980370625e+05, 1.5489835e+06}];

core_depths = [1.38, 1.73]; % in m

%% load lakes

[lakes, exp_volumes] = load_lakes_volumes(lake_defs, core_depths);
num_lakes = size(lakes,1);

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

