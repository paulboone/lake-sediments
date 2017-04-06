function [lakes, exp_volumes] = load_lakes_volumes(lake_defs, core_depths)

  num_lakes = size(lake_defs,1);

  lakes(1,num_lakes) = Lake();
  exp_volumes = zeros(num_lakes);
  for i = 1:num_lakes
    lake_def = lake_defs(i,:);
    lakes(i).load_from_geotiff(lake_def{:});
    exp_volumes(i) = lakes(i).calculate_sediment_volume_from_core(core_depths(i));
  end
end