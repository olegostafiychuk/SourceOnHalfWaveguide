function [a_p_field_1_forw, a_p_field_2_forw, psi_forward_1, psi_forward_2,...
    B_1_forward_1,B_2_forward_1,Cm2_forward_1,  Dm2_forward_1,...
    B_1_forward_2,B_2_forward_2,Cm2_forward_2,  Dm2_forward_2] = coefficientsOf_deltaFunctionWave(typeOfCylinder, q, p, k_0, k, a_0, EE1, GG1, HH1, MU1, EE, MU, m, AE_0, AH_0)

if(abs(real(p))<=1 && imag(q)== 0)
        %%% first type of waves
        if(strcmp(typeOfCylinder, 'Isotropic'))
            psi_forward_1  = psi1_q(k_0, k, a_0, EE1, MU1, EE, MU, m,   p, q);
            psi_backward_1 = psi1_q(k_0, k, a_0, EE1, MU1, EE, MU, m, (-p), q);
        elseif(strcmp(typeOfCylinder, 'Gyrotropic'))
            psi_forward_1  = psi1_q__gyrotropic(k_0, k, a_0, EE1, GG1, HH1, m,   p, q);
            psi_backward_1 = psi1_q__gyrotropic(k_0, k, a_0, EE1, GG1, HH1, m, (-p), q);
        end
        [B_1_forward_1,B_2_forward_1,Cm2_forward_1,  Dm2_forward_1]  = coefficientsOfContinuousSpectrum(typeOfCylinder,k_0, k, a_0, EE1, GG1, HH1, MU1, EE, MU, m,   p, q,   psi_forward_1,  1);
        [B_1_backward_1,B_2_backward_1,Cm2_backward_1, Dm2_backward_1] = coefficientsOfContinuousSpectrum(typeOfCylinder,k_0, k, a_0, EE1, GG1, HH1, MU1, EE, MU, m, (-p), q, psi_backward_1, 1);
        [B_1_back_PandM_1,B_2_back_PandM_1,Cm2_back_PandM_1, Dm2_back_PandM_1] = coefficientsOfContinuousSpectrum(typeOfCylinder,k_0, k, a_0, EE1, -GG1, HH1, MU1, EE, MU, -m, (-p), q, psi_backward_1, 1);

        %%% second type of waves
        if(strcmp(typeOfCylinder, 'Isotropic'))
            psi_forward_2  = psi2_q(k_0, k, a_0, EE1, MU1, EE, MU, m,   p, q);
            psi_backward_2 = psi2_q(k_0, k, a_0, EE1, MU1, EE, MU, m, (-p), q);
        elseif(strcmp(typeOfCylinder, 'Gyrotropic'))
            psi_forward_2  = psi2_q__gyrotropic(k_0, k, a_0, EE1, GG1, HH1, m,   p, q);
            psi_backward_2 = psi2_q__gyrotropic(k_0, k, a_0, EE1, GG1, HH1, m, (-p), q);
        end
        [B_1_forward_2,B_2_forward_2,Cm2_forward_2,  Dm2_forward_2]  = coefficientsOfContinuousSpectrum_forPsi2(typeOfCylinder,k_0, k, a_0, EE1, GG1, HH1, MU1, EE, MU, m,   p, q,   psi_forward_2,  1);
        [B_1_backward_2,B_2_backward_2,Cm2_backward_2, Dm2_backward_2] = coefficientsOfContinuousSpectrum_forPsi2(typeOfCylinder,k_0, k, a_0, EE1, GG1, HH1, MU1, EE, MU, m, (-p), q, psi_backward_2, 1);
        [B_1_back_PandM_2,B_2_back_PandM_2,Cm2_back_PandM_2, Dm2_back_PandM_2] = coefficientsOfContinuousSpectrum_forPsi2(typeOfCylinder,k_0, k, a_0, EE1, -GG1, HH1, MU1, EE, MU, -m, (-p), q, psi_backward_2, 1);

        
        %%%%% forward waves first type of waves %%%%%%%%%%%%%%%%%%%%%%%%%%%
        N_p_1 = - (16 * pi / k_0^2) * (p./ q).* (Cm2_forward_1.^2 + Dm2_forward_1.^2).* psi_forward_1;
        a_p_field_1_forw = - (8 * pi^2 / k_0^3) * (p./ q.^2).*...
                             (Cm2_back_PandM_1.* (psi_forward_1 + 1).* AE_0 +...
                              Dm2_back_PandM_1.* (psi_forward_1 - 1).* AH_0)./ N_p_1;
                                           
        %%%%% forward waves second type of waves %%%%%%%%%%%%%%%%%%%%%%%%%%
        N_p_2 = - (16 * pi / k_0^2) * (p./ q).* (Cm2_forward_2.^2 + Dm2_forward_2.^2).* psi_forward_2;
        a_p_field_2_forw = - (8 * pi^2 / k_0^3) * (p./ q.^2).*...
                             (Cm2_back_PandM_2.* (psi_forward_2 + 1).* AE_0 +...
                              Dm2_back_PandM_2.* (psi_forward_2 - 1).* AH_0)./ N_p_2;
else
    %%% discret mode
        %%% first type of waves
        psi_forward_1  = 0;
        psi_backward_1 = 0;

        [B_1_forward_1,B_2_forward_1,Cm2_forward_1,  Dm2_forward_1]  = coefficientsOfContinuousSpectrum(typeOfCylinder,k_0, k, a_0, EE1, GG1, HH1, MU1, EE, MU, m,   p, q,   psi_forward_1,  1);
        [B_1_backward_1,B_2_backward_1,Cm2_backward_1, Dm2_backward_1] = coefficientsOfContinuousSpectrum(typeOfCylinder,k_0, k, a_0, EE1, GG1, HH1, MU1, EE, MU, m, (-p), q, psi_backward_1, 1);
        [B_1_back_PandM_1,B_2_back_PandM_1,Cm2_back_PandM_1, Dm2_back_PandM_1] = coefficientsOfContinuousSpectrum(typeOfCylinder,k_0, k, a_0, EE1, -GG1, HH1, MU1, EE, MU, -m, (-p), q, psi_backward_1, 1);

        %%% second type of waves
        psi_forward_2  = 0;
        psi_backward_2 = 0;

        [B_1_forward_2,B_2_forward_2,Cm2_forward_2,  Dm2_forward_2]  = coefficientsOfContinuousSpectrum_forPsi2(typeOfCylinder,k_0, k, a_0, EE1, GG1, HH1, MU1, EE, MU, m,   p, q,   psi_forward_2,  1);
        [B_1_backward_2,B_2_backward_2,Cm2_backward_2, Dm2_backward_2] = coefficientsOfContinuousSpectrum_forPsi2(typeOfCylinder,k_0, k, a_0, EE1, GG1, HH1, MU1, EE, MU, m, (-p), q, psi_backward_2, 1);
        [B_1_back_PandM_2,B_2_back_PandM_2,Cm2_back_PandM_2, Dm2_back_PandM_2] = coefficientsOfContinuousSpectrum_forPsi2(typeOfCylinder,k_0, k, a_0, EE1, -GG1, HH1, MU1, EE, MU, -m, (-p), q, psi_backward_2, 1);

        
        %%%%% forward waves first type of waves %%%%%%%%%%%%%%%%%%%%%%%%%%%
        a_p_field_1_forw = 1;
                                           
        %%%%% forward waves second type of waves %%%%%%%%%%%%%%%%%%%%%%%%%%
        a_p_field_2_forw = 0;
                          
end
                          
                                                                                 
