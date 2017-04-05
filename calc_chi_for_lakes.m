function [chi2, chi2all] = calc_chi_for_lakes(lakes, exp_volumes, ttlem_params, k, d)
  ttlem_params.D = d;
  ttlem_params.Kw = k;

  ttlem_params = ttlemset(ttlem_params);

  chi2 = zeros(1,length(lakes));
  for l=1:length(lakes)
    v_mod = lakes(l).calculate_sediment_volume_via_model(ttlem_params);
    chi2(l) = ((v_mod - exp_volumes(l))^2)/(exp_volumes(l)^2);
  end

  chi2all = squeeze(sum(chi2));

  close all;
end