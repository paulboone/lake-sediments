
close all, clear all, clc
global lakes exp_volumes m

%% configuration
matfilename = 'all_k0.000625_0.00125_d0.05625_0.09375.mat';

% run optimizations from the {num_starting_bins} bins with the lowest chi^2
% values.
num_starting_bins = 3;

% or set the above to NaN and manually set the start k/d vals:
% search_kd = [0.0009063,0.0775];

% can set parameters to fminsearch here via optim set
search_options = optimset();

% search_options = optimset('MaxFunEvals',1); % for debugging



%% run 

load(matfilename, 'm', 'lakechi2');

% find bins with lowest chi2 to start in:
if ~isnan(num_starting_bins)
  [~,sort_indices] = sort(lakechi2(:),'ascend');
  search_indices = sort_indices(1:num_starting_bins);

  search_kd = zeros(num_starting_bins,2);
  for i = 1:num_starting_bins
    [vi,ki] = ind2sub(size(lakechi2),search_indices(i));
    search_kd(i,:) = [m.k_vals(ki),m.d_vals(vi)];
  end
end



[lakes, exp_volumes] = load_lakes_volumes(m.lake_defs, m.core_depths);
num_lakes = size(lakes,1);

mins = search_kd * 0;
chi2 = nan(1,2);
for i = 1:size(search_kd,1)
  [mins(i,1:2), chi2(i)] = fminsearch(@sediment_volume_kd, search_kd(i,:),search_options);
end

%% plot
plot_chi2('Lake', lakechi2, m.k_vals, m.d_vals);
add_optimization_points_to_chi2_plot(mins, chi2);
export_fig(['chi2_', matfilename,  'optim.png']);


function chi2 = sediment_volume_kd(kd)
  global lakes exp_volumes m

  k = kd(1);
  d = kd(2);

  if k < 0 || d < 0
    chi2 = 1000;
  else
    [~, chi2] = calc_chi_for_lakes(lakes, exp_volumes, m.ttlem_params, k, d);
  end

  disp([k, d, chi2]);

  close all;
end
