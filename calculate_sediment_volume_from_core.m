function sediment_volume = calculate_sediment_volume_from_core(lake, sediment_depth)
  sediment_volume = sediment_depth * lake.num_lake_cells * lake.cell_area;
end
