%%%%% расчёт коэффициента для волн непрерывного спектра в праямом
%%%%% направлении (в сторону среды 2)
%%%%% второе приближение в нахождении коэффициента aplus
function aplus2 = a_sma_of_continuousWaves_alpha1(typeOfCylinder, q, p_0, p, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, z, I_f, I_z, d)

switch(typeOfCylinder)
    case 'Isotropic'
        psi  = psi1_q(k_0, k_0, a_0, EE1, MU1, EE, MU, m, p, q);
        [B_1, B_2, Cm2, Dm2]  = coefficientsOfContinuousSpectrum(typeOfCylinder, k_0, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, p, q, psi,  1);
        %%%%%% norm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        N_p = Norm_of_continuousWaves(q, p, k_0, c, m, psi, Cm2, Dm2);
        
        psi  = psi1_q(k_0, k_0, a_0, EE1, MU1, EE, MU, -m, -p, q);
        [B_1, B_2, Cm2, Dm2]  = coefficientsOfContinuousSpectrum(typeOfCylinder, k_0, k_0, a_0, EE1, -GG1, HH1, MU1, EE, MU, -m, -p, q, psi,  1);
        
        Ez  = field_continuesWaves_discreteRepr('Ez', typeOfCylinder, a_0, q, -p, k_0, a_0, EE1, -GG1, HH1, MU1, -m, z,...
          1, psi, B_1,B_2,Cm2,  Dm2);
        Ephi= field_continuesWaves_discreteRepr('Ephi', typeOfCylinder, a_0, q, -p, k_0, a_0, EE1,-GG1, HH1, MU1,-m, z,...
          1, psi, B_1,B_2,Cm2,  Dm2);
        
        aplus2 = (2 * pi * a_0./ N_p).* (I_z.* Ez + I_f.* Ephi).* 2.* sin(k_0 * (p - p_0) * d)./ (k_0 * (p - p_0));
    
         
         
    case 'Gyrotropic'
        psi   = psi1_q__gyrotropic(k_0, k_0, a_0, EE1,  GG1, HH1,  m, p, q);
        [B_1, B_2, Cm2, Dm2]  = coefficientsOfContinuousSpectrum(typeOfCylinder, k_0, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, p, q, psi,  1);
        %%%%%% norm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        N_p = Norm_of_continuousWaves(q, p, k_0, c, m, psi, Cm2, Dm2);
        
        psi   = psi1_q__gyrotropic(k_0, k_0, a_0, EE1,  -GG1, HH1,  -m, -p, q);
        [B_1, B_2, Cm2, Dm2]  = coefficientsOfContinuousSpectrum(typeOfCylinder, k_0, k_0, a_0, EE1, -GG1, HH1, MU1, EE, MU, -m, -p,  q, psi,  1);
        
        Ez  = field_continuesWaves_discreteRepr('Ez', typeOfCylinder, a_0 , q, -p, k_0, a_0, EE1, -GG1, HH1, MU1, -m, z,...
          1, psi, B_1,B_2,Cm2,  Dm2);
        Ephi= field_continuesWaves_discreteRepr('Ephi', typeOfCylinder, a_0, q, -p, k_0, a_0, EE1, -GG1, HH1, MU1, -m, z,...
          1, psi, B_1,B_2,Cm2,  Dm2);

       aplus2 = (2 * pi * a_0./ N_p).* (I_z.* Ez + I_f.* Ephi).* 2.* sin(k_0 * (p - p_0) * d)./ (k_0 * (p - p_0));
end


    
    
    
    
    
    
    
    
    
