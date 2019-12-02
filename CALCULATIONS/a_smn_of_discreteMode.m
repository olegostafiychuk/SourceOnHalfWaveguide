%%%%% расчёт коэффициента для волн непрерывного спектра в праямом
%%%%% направлении (в сторону среды 2)
%%%%% второе приближение в нахождении коэффициента aplus
function aplus2 = a_smn_of_discreteMode(typeOfCylinder, q, p_0, p, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, j_f, j_z, d)

  switch(typeOfCylinder)
      case 'Isotropic'  
          N_n = Norm_of_descreteMode(typeOfCylinder, q, p, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m);

          Ez  = field_of_discreteMode('Ez',   typeOfCylinder, a_0, q, p, 1, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, 0);
          Ephi= field_of_discreteMode('Ephi', typeOfCylinder, a_0, q, p, 1, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, 0);
          aplus2 = (2 * pi * a_0./ N_n).* (j_z.* Ez + j_f.* Ephi).* 2.* sin(k_0 * (p - p_0) * d)./ (k_0 * (p - p_0));
          
      case 'Gyrotropic'
          N_n = Norm_of_descreteMode(typeOfCylinder, q, p, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m);
          
          Ez  = field_of_discreteMode('Ez',   typeOfCylinder, a_0, q, -p, 1, k_0, a_0, EE1, -GG1, HH1, MU1, EE, MU, -m, 0);
          Ephi= field_of_discreteMode('Ephi', typeOfCylinder, a_0, q, -p, 1, k_0, a_0, EE1, -GG1, HH1, MU1, EE, MU, -m, 0);
          
          aplus2 = (2 * pi * a_0./ N_n).* (j_z.* Ez + j_f.* Ephi).* 2.* sin(k_0 * (p - p_0) * d)./ (k_0 * (p - p_0));
  end
