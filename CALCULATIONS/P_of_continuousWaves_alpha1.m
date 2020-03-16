%%%%% расчёт коэффициента для волн непрерывного спектра в праямом
%%%%% направлении (в сторону среды 2)
%%%%% второе приближение в нахождении коэффициента aplus
function P = P_of_continuousWaves_alpha1(typeOfCylinder, q, p_0, p, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, z, I_f, I_z, d, aplus2)

switch(typeOfCylinder)
    case 'Isotropic'
        psi  = psi1_q(k_0, k_0, a_0, EE1, MU1, EE, MU, m, p, q);
        [B_1, B_2, Cm2, Dm2]  = coefficientsOfContinuousSpectrum(typeOfCylinder, k_0, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, p, q, psi,  1);
        %%%%%% norm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        N_p = Norm_of_continuousWaves(q, p, k_0, c, m, psi, Cm2, Dm2);
      
        %aplus2 = a_sma_of_continuousWaves_alpha1(typeOfCylinder, q, p_0, p, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, z, I_f, I_z, d);
        
        P = 1/4 * abs(aplus2).^2.* N_p;
         
    case 'Gyrotropic'
        psi   = psi1_q__gyrotropic(k_0, k_0, a_0, EE1,  GG1, HH1,  m, p, q);
        [B_1, B_2, Cm2, Dm2]  = coefficientsOfContinuousSpectrum(typeOfCylinder, k_0, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, p, q, psi,  1);
        %%%%%% norm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        N_p = Norm_of_continuousWaves(q, p, k_0, c, m, psi, Cm2, Dm2);       
%         aplus2 = a_sma_of_continuousWaves_alpha1(typeOfCylinder, q, p_0, p, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, z, I_f, I_z, d);        
        %P = 1/4 * real(abs(aplus2).^2.* N_p); 
        
%         P_sm1 = c / k_0^2 * p./ q.* (abs(Cm2).^2 + abs(Dm2).^2).* psi;
        P_sm1 = c /(2*k_0^2) * p./ q.* (abs(Cm2).^2 + abs(Dm2).^2).* (1 + abs(psi).^2);
        P = real(abs(aplus2).^2.* P_sm1);
end


    
    
    
    
    
    
    
    
    
