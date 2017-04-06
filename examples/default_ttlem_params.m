
function ttlem_params = default_ttlem_params()
  ttlem_params.TimeSpan = 6000;
  ttlem_params.TimeStep = 6000;
  ttlem_params.m = 0.5;
  ttlem_params.n = 1;
  ttlem_params.DrainDir = 'variable';
  ttlem_params.AreaThresh = 2e5;     % channel contributing area threshold
  ttlem_params.Sc = 1;
  ttlem_params.Sc_unit = 'tangent';
  ttlem_params.ploteach=inf;
  ttlem_params.saveeach=1;
end