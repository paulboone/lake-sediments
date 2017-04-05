function plotsave_chi2(lake_name, lakechi2, k_vals, d_vals)
  figure
  imagesc(k_vals, d_vals, log10(lakechi2))
  title(['Chi^2 graph of varying K,D values (', lake_name, ')'])
  xlabel('K [m^2yr^{-1}]');
  ylabel('D [m^2yr^{-1}]');
  xticks(k_vals);
  yticks(d_vals);
  caxis([-6, 0])
  colorbar
  shg
  file_basename = [lake_name,'_k', num2str(k_vals(1)), '_', num2str(k_vals(end)), '_d', num2str(d_vals(1)), '_', num2str(d_vals(end))];
  save([file_basename, '.mat'], 'k_vals', 'd_vals', 'lakechi2');
  export_fig(['chi2_', file_basename,  '.png']);
end
