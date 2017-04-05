
load('all_k0.000875_0.001125_d0.07_0.08.mat', 'd_vals', 'k_vals', 'lakechi2');

plot_chi2('Lake', lakechi2, k_vals, d_vals);

optx = [0.0009, 0.00095];
opty = [0.075, 0.077];

add_optimization_points_to_chi2_plot(optx, opty)

function add_optimization_points_to_chi2_plot(optx, opty)
  hold on;
  
  plot(optx, opty, 'o', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');

  for i = 1:length(optx)
    text_label = char(strcat('  - (', string(optx(i)), ', ', string(opty(i)),')'));
    text(optx(i), opty(i),text_label);
  end

  hold off;
end
