function plotsave_chi2(lake_name, lakechi2, k_vals, d_vals)
  plot_chi2(lake_name, lakechi2, k_vals, d_vals);
  shg
  file_basename = [lake_name,'_k', num2str(k_vals(1)), '_', num2str(k_vals(end)), '_d', num2str(d_vals(1)), '_', num2str(d_vals(end))];
  save([file_basename, '.mat'], 'k_vals', 'd_vals', 'lakechi2');
  export_fig(['chi2_', file_basename,  '.png']);
end
