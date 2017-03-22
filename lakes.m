
% add topotoolbox paths - double check what is needed!
addpath ~/workspace/topotoolbox/
addpath ~/workspace/topotoolbox/utilities
addpath ~/workspace/topotoolbox/topoapp
addpath ~/workspace/topotoolbox/DEMdata
addpath ~/workspace/topotoolbox/ttlem
addpath ~/workspace/topotoolbox/ttlem/LateralDisplacement
addpath ~/workspace/export_fig
addpath ~/workspace/freezeColors


close all, clear all, clc



%% TTLEM vars
ttlem_params.TimeSpan = 6000;
ttlem_params.TimeStep = 6000;

% Hillslope processes, parameters values
ttlem_params.m = 0.5;
ttlem_params.n = 1;

% Fixed or variable drainage network through time
ttlem_params.DrainDir = 'variable';

ttlem_params.AreaThresh = 2e5;     % channel contributing area threshold, m�

% Threshold slopes
ttlem_params.Sc = 1;
ttlem_params.Sc_unit = 'tangent';

% Output
ttlem_params.ploteach=inf;
ttlem_params.saveeach=1;

%% run lake

ttlem_params.D = 0.1;
ttlem_params.Kw = 0.05;

project_dir = '/Users/pboone/Dropbox/Projects/Classes/GEOL 2049/project-lake-sediment-analysis/';
lost_lake = init_lake([project_dir, 'clipped_lost.tif'], 3.523578530974025e+05, 1.612583134776515e+06);
jan_lake = init_lake([project_dir, 'clipped_jan1.tif'], 4.980370625e+05, 1.5489835e+06);


v_exp_jan = calculate_sediment_volume_from_core(jan_lake, 1.73);
v_exp_lost = calculate_sediment_volume_from_core(lost_lake, 1.38);
% v_mod = calculate_sediment_volume_via_model(lost_lake, ttlem_params);


% v_mod = calculate_sediment_volume_via_model(jan_lake, ttlem_params);

% d_vals = linspace(0, 0.3, 10);    %0:0.05:2.5*10
% k_vals = linspace(0, 0.01, 10);  %3e-5:0.0005:3e-2
d_vals = linspace(0, 0.5, 20);    %0:0.05:2.5*10
k_vals = linspace(0, 0.025, 20);  %3e-5:0.0005:3e-2


% D=linspace(0, 0.5,20);      %0:0.05:2.5*10
% K=linspace(0, 0.025,20);    %3e-5:0.0005:3e-2

% [Md, Mk]=meshgrid(D,K);

for i=1:length(d_vals);
  for j=1:length(k_vals)
    ttlem_params.D = d_vals(i);
    ttlem_params.Kw = k_vals(j);

    ttlem_params = ttlemset(ttlem_params);
    v_mod = calculate_sediment_volume_via_model(jan_lake, ttlem_params);
    chi2(i,j) = ((v_mod - v_exp_jan)^2)/(v_exp_jan^2);

    % v_mod = calculate_sediment_volume_via_model(lost_lake, ttlem_params);
    % chi2(i,j) = chi2(i,j) + ((v_mod - v_exp_lost)^2)/(v_exp_lost^2);

    close all;

    disp(ttlem_params);
    disp([i, j, v_mod, chi2(i,j)]);
  end
  disp(i);
end

imagesc(k_vals, d_vals, log10(chi2)), shg
xlabel('K [m^{1-2m}yr^{-1}]');
ylabel('D [m^2yr^{-1}]');



%% sediment_volume(sediment_depth)
%% calculate sediment volume from lake size
function sediment_volume = calculate_sediment_volume_from_core(lake, sediment_depth)
  sediment_volume = sediment_depth * lake.num_lake_cells * lake.cell_area;
end
