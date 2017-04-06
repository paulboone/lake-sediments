function plotsave_chi2(lake_name, lakechi2, m)
  plot_chi2(lake_name, lakechi2, m.k_vals, m.d_vals);
  shg
  file_basename = [lake_name,'_k', num2str(m.k_vals(1)), '_', num2str(m.k_vals(end)), '_d', num2str(m.d_vals(1)), '_', num2str(m.d_vals(end))];
  save([file_basename, '.mat'], 'lakechi2','m');
  export_fig(['chi2_', file_basename,  '.png']);
end
