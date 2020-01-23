%%%%% расчёт коэффициента для волн непрерывного спектра в праямом
%%%%% направлении (в сторону среды 2)
%%%%% второе приближение в нахождении коэффициента aplus
function aplus2 = a_p_field__ForLGBeam(typeOfCylinder, q, q_0, p, k_0, k, a_0, EE1, GG1, HH1, MU1, EE2, MU2, EE, MU, m, z, psi, a_b, LG_AE_0, LG_AH_0, upperBound)

switch(typeOfCylinder)
    case 'Gyrotropic'
      
      
        
 %%%%%%%%%%%%%%%%%%%%% outer integral %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [q1, q2, n1, n2, alp1, alp2, bet1, bet2] = term_of_gyrotropic_waveguide(EE1, -GG1, HH1, -p);
        [B_1, B_2, Cm2, Dm2]  = coefficientsOfContinuousSpectrum(typeOfCylinder, k_0, k, a_0, EE1, -GG1, HH1, MU1, EE, MU, -m, -p, q, psi,  1);
        %%%%%% norm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        N_p = -(-1)^m * (8 / k_0) * (p./ q).* (Cm2.^2 + Dm2.^2).* psi;
        
for iq = 1:size(q,2)
%     if(real(q(iq)) < 100)
%         if(q(iq) > 10)
%         S_plus(iq)  = quadgk(@(r) r.* L_core_gyrWaveguide_outer(r, q(iq), -p(iq), -m, k_0, a_0, m, a_b, LG_AE_0, LG_AH_0, Cm2(iq), Dm2(iq), psi(iq)), a_0, upperBound * 100/q(iq));
%     else
        S_plus(iq)  = quadgk(@(r) r.* L_core_gyrWaveguide_outer(r, q(iq), -p(iq), -m, k_0, a_0, m, a_b, LG_AE_0, LG_AH_0, Cm2(iq), Dm2(iq), psi(iq)), a_0, upperBound);
%     end
%         q(iq)
    
        S_minus(iq) = quadgk(@(r) r.* L_core_gyrWaveguide_inner(r, q(iq), -p(iq), -m, k_0, a_0, m, a_b, LG_AE_0, LG_AH_0,...
                             B_1(iq), B_2(iq), q1(iq), q2(iq), n1(iq), n2(iq), alp1(iq), alp2(iq), bet1(iq), bet2(iq), EE1, -GG1, HH1), 0, a_0);
%     else
%         S_plus(iq) = 0; S_minus(iq) = 0;
%     end
        
end
                                                                               
        
       aplus2 = (S_plus + S_minus)./ N_p;
end


 
