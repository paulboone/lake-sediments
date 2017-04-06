function add_optimization_points_to_chi2_plot(mins, chi2)
  hold on;

  plot(mins(:,1), mins(:,2), 'o', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');

  legend_labels = strings(size(mins,1));
  for i = 1:size(mins,1)
    legend_labels(i) = char(strcat(num2str(i),') \chi= ', num2str(chi2(i),'%4.3e'),' @ (', num2str(mins(i,1)), ', ', num2str(mins(i,2)),')'));
    text_label = strcat(' - (', num2str(i), ')');
    text(mins(i,1), mins(i,2),text_label);
  end
  
  legend(strjoin(legend_labels, '\n'),'Location','southoutside');

  hold off;
end