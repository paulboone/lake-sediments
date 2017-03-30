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

global lake ttlem_params v_exp

lake = Lake();
lake.load_from_geotiff('clipped_jan.tif', NaN, NaN);
v_exp = lake.calculate_sediment_volume_from_core(1.73);

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

xmin = fminsearch(@sediment_volume_kd, [0.0011, 0.06]);
disp(xmin);

function chi2 = sediment_volume_kd(kd)
  global lake ttlem_params v_exp
  
  k = kd(1);
  d = kd(2);
  
  if k < 0 || d < 0
    chi2 = 1000;
  else
    ttlem_params.D = d;
    ttlem_params.Kw = k;
    ttlem_params = ttlemset(ttlem_params);

    v_mod = lake.calculate_sediment_volume_via_model(ttlem_params);
    chi2 = ((v_mod - v_exp)^2)/(v_exp^2);
  end
  
  disp([k, d, chi2]);
    
  close all;
end

