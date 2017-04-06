function plot_chi2(lake_name, lakechi2, k_vals, d_vals)
  figure('Position', [0, 0, 800, 800])
  imagesc(k_vals, d_vals, log10(lakechi2))
  title(['\chi^2 graph of varying K,D values (', lake_name, ')'])
  xlabel('K [m^2yr^{-1}]');
  ylabel('D [m^2yr^{-1}]');
  xticks(k_vals);
  yticks(d_vals);
  caxis([-6, 0])
  colorbar
  
  shg
end
