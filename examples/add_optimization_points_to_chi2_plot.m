function add_optimization_points_to_chi2_plot(mins, chi2)
  hold on;

  plot(mins(:,1), mins(:,2), 'o', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');

  for i = 1:size(mins,1)
    text_label = char(strcat(' - \chi= ', string(chi2(i)),' @ (', string(mins(i,1)), ', ', string(mins(i,2)),')'));
    text(mins(i,1), mins(i,2),text_label);
  end

  hold off;
end