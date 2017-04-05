function plot_chi2(lake_name, lakechi2, k_vals, d_vals)
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
end
