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

%% setup params
m.ttlem_params = default_ttlem_params();

m.lake_defs = [{'geotiffs/lost.tif', 3.525e+05, 1.6125e+06}; ...
               {'geotiffs/jan.tif', 4.98e5, 1.549e6}];
m.core_depths = [1.38, 1.73]; % in m

m.k_vals = linspace(0.0, 0.01, 9);
m.d_vals = linspace(0.0, 0.3, 9);

%% run lakes
[lakes, exp_volumes] = load_lakes_volumes(m.lake_defs, m.core_depths);
num_lakes = length(lakes);

chi2 = NaN(num_lakes, length(m.d_vals), length(m.k_vals));
chi2all = NaN(length(m.d_vals), length(m.k_vals));
for i=1:length(m.d_vals)
  disp(['Progress: ' num2str(i-1) '/' num2str(length(m.d_vals))]);
  for j=1:length(m.k_vals)
    [chi2(:,i,j), chi2all(i,j)] = calc_chi_for_lakes(lakes, exp_volumes, m.ttlem_params, m.k_vals(j), m.d_vals(i));
  end
end

%% plot
for i=1:num_lakes
  lake = lakes(i);
  m_lakeonly = m;
  m_lakeonly.lake_defs = m.lake_defs(i,:);
  m_lakeonly.core_depths = m.core_depths(i);
  
  lakechi2 = squeeze(chi2(i,:,:));
  plotsave_chi2(lake.lake_name, lakechi2, m_lakeonly);
end

plotsave_chi2('all', chi2all, m);
