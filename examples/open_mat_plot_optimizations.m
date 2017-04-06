
global lakes ttlem_params exp_volumes

matfilename = 'all_k0.000875_0.001125_d0.07_0.08.mat';

lake_defs = [{'clipped_lost.tif', 3.523578530974025e+05, 1.612583134776515e+06}; ...
             {'clipped_jan.tif', 4.980370625e+05, 1.5489835e+06}];

core_depths = [1.38, 1.73]; % in m
search_kd = [0.0009063,0.0775]; %; 0.001, 0.0738
ttlem_params = default_ttlem_params();

search_options = optimset();
%% run optimizations

[lakes, exp_volumes] = load_lakes_volumes(lake_defs, core_depths);
num_lakes = size(lakes,1);

mins = search_kd * 0;
chi2 = nan(1,2);
for i = 1:size(search_kd,1)
  [mins(i,1:2), chi2(i)] = fminsearch(@sediment_volume_kd, search_kd(i,:),search_options);
end


%% plot
load(matfilename', 'd_vals', 'k_vals', 'lakechi2');
plot_chi2('Lake', lakechi2, k_vals, d_vals);
add_optimization_points_to_chi2_plot(mins, chi2)
export_fig(['chi2_', matfilename,  'optim.png']);


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
