function [a_p_field_Ewave_forw, a_p_field_Hwave_forw] = coefficientsOf_deltaFunctionWaveOfFreeSpace(typeOfCylinder, q, p, k_0, k, a_0, EE1, GG1, HH1, MU1, EE, MU, m,...
    a_p_field_1_forw, a_p_field_2_forw, psi_forward_1, psi_forward_2,...
    Cm2_forward_1,  Dm2_forward_1,...
    Cm2_forward_2,  Dm2_forward_2)

%         %%% first type of waves
%         if(strcmp(typeOfCylinder, 'Isotropic'))
%             psi_forward_1  = psi1_q(k_0, k, a_0, EE1, MU1, EE, MU, m,   p);
%             psi_backward_1 = psi1_q(k_0, k, a_0, EE1, MU1, EE, MU, m, (-p));
%         elseif(strcmp(typeOfCylinder, 'Gyrotropic'))
%             psi_forward_1  = psi1_q__gyrotropic(k_0, k, a_0, EE1, GG1, HH1, m,   p);
%             psi_backward_1 = psi1_q__gyrotropic(k_0, k, a_0, EE1, GG1, HH1, m, (-p));
%         end
%         [B_1_forward_1,B_2_forward_1,Cm2_forward_1,  Dm2_forward_1]  = coefficientsOfContinuousSpectrum(typeOfCylinder,k_0, k, a_0, EE1, GG1, HH1, MU1, EE, MU, m,   p,   psi_forward_1,  1);
%         [B_1_backward_1,B_2_backward_1,Cm2_backward_1, Dm2_backward_1] = coefficientsOfContinuousSpectrum(typeOfCylinder,k_0, k, a_0, EE1, GG1, HH1, MU1, EE, MU, m, (-p), psi_backward_1, 1);
%         [B_1_back_PandM_1,B_2_back_PandM_1,Cm2_back_PandM_1, Dm2_back_PandM_1] = coefficientsOfContinuousSpectrum(typeOfCylinder,k_0, k, a_0, EE1, -GG1, HH1, MU1, EE, MU, -m, (-p), psi_backward_1, 1);

%         %%% second type of waves
%         if(strcmp(typeOfCylinder, 'Isotropic'))
%             psi_forward_2  = psi2_q(k_0, k, a_0, EE1, MU1, EE, MU, m,   p);
%             psi_backward_2 = psi2_q(k_0, k, a_0, EE1, MU1, EE, MU, m, (-p));
%         elseif(strcmp(typeOfCylinder, 'Gyrotropic'))
%             psi_forward_2  = psi2_q__gyrotropic(k_0, k, a_0, EE1, GG1, HH1, m,   p);
%             psi_backward_2 = psi2_q__gyrotropic(k_0, k, a_0, EE1, GG1, HH1, m, (-p));
%         end
%         [B_1_forward_2,B_2_forward_2,Cm2_forward_2,  Dm2_forward_2]  = coefficientsOfContinuousSpectrum(typeOfCylinder,k_0, k, a_0, EE1, GG1, HH1, MU1, EE, MU, m,   p,   psi_forward_2,  1);
%         [B_1_backward_2,B_2_backward_2,Cm2_backward_2, Dm2_backward_2] = coefficientsOfContinuousSpectrum(typeOfCylinder,k_0, k, a_0, EE1, GG1, HH1, MU1, EE, MU, m, (-p), psi_backward_2, 1);
%         [B_1_back_PandM_2,B_2_back_PandM_2,Cm2_back_PandM_2, Dm2_back_PandM_2] = coefficientsOfContinuousSpectrum(typeOfCylinder,k_0, k, a_0, EE1, -GG1, HH1, MU1, EE, MU, -m, (-p), psi_backward_2, 1);

        if(p<=1)
        %%%%% forward waves first type of waves %%%%%%%%%%%%%%%%%%%%%%%%%%%
        N_p_1 = - (4 * pi / k_0^2) * (p./ q);
        a_p_field_Ewave_forw = - (4 * pi / k_0^2) * (p./ q).*...
                                              (Cm2_forward_1.* (psi_forward_1 + 1).* a_p_field_1_forw +...
                                               Cm2_forward_2.* (psi_forward_2 + 1).* a_p_field_2_forw)./ N_p_1;
                                           
        %%%%% forward waves second type of waves %%%%%%%%%%%%%%%%%%%%%%%%%%
        N_p_2 =   (4 * pi / k_0^2) * (p./ q);
        a_p_field_Hwave_forw = - (4 * pi / k_0^2) * (p./ q).*...
                                              (Dm2_forward_1.* (psi_forward_1 - 1).* a_p_field_1_forw +...
                                               Dm2_forward_2.* (psi_forward_2 - 1).* a_p_field_2_forw)./ N_p_2;
        else
        %%%%% forward waves first type of waves %%%%%%%%%%%%%%%%%%%%%%%%%%%
        a_p_field_Ewave_forw = 0;
                                           
        %%%%% forward waves second type of waves %%%%%%%%%%%%%%%%%%%%%%%%%%
        a_p_field_Hwave_forw = 0;
        end
                                                                                 
