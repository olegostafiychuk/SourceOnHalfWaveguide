%%%%% расчёт коэффициента для волн непрерывного спектра в прямом
%%%%% направлении (в сторону среды 2)
%%%%% второе приближение в нахождении коэффициента aplus
function [a_p_field_1_forw, a_p_field_1_back, a_p_field_2_forw, a_p_field_2_back,...
          psi_forward_1, psi_backward_1, psi_forward_1_transp, psi_backward_1_transp,...
          psi_forward_2, psi_backward_2, psi_forward_2_transp, psi_backward_2_transp,...
          B_1_forward_1,B_2_forward_1,Cm2_forward_1,  Dm2_forward_1,...
          B_1_backward_1,B_2_backward_1,Cm2_backward_1, Dm2_backward_1,...
          B_1_forward_2,B_2_forward_2,Cm2_forward_2,  Dm2_forward_2,...
          B_1_backward_2,B_2_backward_2,Cm2_backward_2, Dm2_backward_2] =...
          coeffsOfFieldAndExcitCoeffsOfLGBeam_of_Continues_discreteRepr(typeOfCylinder, q, q_0, k_0, k, a_0, EE1, GG1, HH1, MU1, EE2, MU2, EE, MU, m, z, a_b, AE_0, AH_0)

  p = sqrt(1-q.^2);
  p = real(p) - 1i * abs(imag(p));
  
  p_0 = sqrt(1-q_0.^2);
  p_0 = real(p_0) - 1i * abs(imag(p_0));
  
  m0 = m;
  
switch(typeOfCylinder)
    case 'Gyrotropic'      
        %%% вычисляем первый тип волн непрерывного спектра
       
        %%% first type of waves
        psi_forward_1   = psi1_q__gyrotropic(k_0, k, a_0, EE1, GG1, HH1, m,   p, q);
        psi_backward_1  = psi1_q__gyrotropic(k_0, k, a_0, EE1, GG1, HH1, m,  -p, q);
        psi_forward_1_transp   = psi1_q__gyrotropic(k_0, k, a_0, EE1, -GG1, HH1, -m,   p, q);
        psi_backward_1_transp  = psi1_q__gyrotropic(k_0, k, a_0, EE1, -GG1, HH1, -m,  -p, q);
%         psi_forward_1_transp   = psi_backward_1;
%         psi_backward_1_transp  = psi_forward_1;
        [B_1_forward_1,B_2_forward_1,Cm2_forward_1,  Dm2_forward_1]  = coefficientsOfContinuousSpectrum(typeOfCylinder,k_0, k, a_0, EE1, GG1, HH1, MU1, EE, MU, m,   p, q,   psi_forward_1, 1);
        [B_1_backward_1,B_2_backward_1,Cm2_backward_1, Dm2_backward_1] = coefficientsOfContinuousSpectrum(typeOfCylinder,k_0, k, a_0, EE1, GG1, HH1, MU1, EE, MU, m, (-p), q, psi_backward_1, 1);
%         [B_1_back_PandM_1,B_2_back_PandM_1,Cm2_back_PandM_1, Dm2_back_PandM_1] = coefficientsOfContinuousSpectrum(typeOfCylinder,k_0, k, a_0, EE1, MU1, EE, MU, -m, (-p), q, psi_forward_1, 1);
        
        %%% second type of waves
        psi_forward_2  = psi2_q__gyrotropic(k_0, k, a_0, EE1, GG1, HH1, m,   p, q);
        psi_backward_2  = psi2_q__gyrotropic(k_0, k, a_0, EE1, GG1, HH1, m,  -p, q);
        psi_forward_2_transp   = psi2_q__gyrotropic(k_0, k, a_0, EE1, -GG1, HH1, -m,   p, q);
        psi_backward_2_transp  = psi2_q__gyrotropic(k_0, k, a_0, EE1, -GG1, HH1, -m,  -p, q);
%         psi_forward_2_transp   = psi_backward_2;
%         psi_backward_2_transp  = psi_forward_2;
        [B_1_forward_2,B_2_forward_2,Cm2_forward_2,  Dm2_forward_2]  =   coefficientsOfContinuousSpectrum_forPsi2(typeOfCylinder,k_0, k, a_0, EE1, GG1, HH1, MU1, EE, MU, m,   p,  q,  psi_forward_2, 1);
        [B_1_backward_2,B_2_backward_2,Cm2_backward_2, Dm2_backward_2] = coefficientsOfContinuousSpectrum_forPsi2(typeOfCylinder,k_0, k, a_0, EE1, GG1, HH1, MU1, EE, MU, m, (-p), q, psi_backward_2, 1);
% %         [B_1_back_PandM_2,B_2_back_PandM_2,Cm2_back_PandM_2, Dm2_back_PandM_2] = coefficientsOfContinuousSpectrum_forPsi2(typeOfCylinder,k_0, k, a_0, EE1, MU1, EE, MU, -m, (-p), q, psi_forward_2, 1);


upperBound = 5*a_b;
        %%%%% forward waves first type of waves %%%%%%%%%%%%%%%%%%%%%%%%%%%
        a_p_field_1_forw = a_p_field__ForLGBeam(typeOfCylinder, q, q_0,  p, k_0, k, a_0, EE1, GG1, HH1,...
            MU1, EE2, MU2, EE, MU, m, z,psi_backward_1_transp, a_b, AE_0, AH_0, upperBound);
        %%%%% backward waves first type of waves %%%%%%%%%%%%%%%%%%%%%%%%%%%
        a_p_field_1_back = a_p_field__ForLGBeam(typeOfCylinder, q, q_0, -p, k_0, k, a_0, EE1, GG1, HH1,...
            MU1, EE2, MU2, EE, MU, m, z,psi_forward_1_transp, a_b, AE_0, AH_0, upperBound);

        
        %%%%% forward waves second type of waves %%%%%%%%%%%%%%%%%%%%%%%%%%%
        a_p_field_2_forw = a_p_field__ForLGBeam_forPsi2(typeOfCylinder, q, q_0,  p, k_0, k, a_0, EE1, GG1, HH1,...
            MU1, EE2, MU2, EE, MU, m, z, psi_backward_2_transp, a_b, AE_0, AH_0, upperBound);
        %%%%% backward waves second type of waves %%%%%%%%%%%%%%%%%%%%%%%%%%%
        a_p_field_2_back = a_p_field__ForLGBeam_forPsi2(typeOfCylinder, q, q_0, -p, k_0, k, a_0, EE1, GG1, HH1,...
            MU1, EE2, MU2, EE, MU, m, z, psi_forward_2_transp, a_b, AE_0, AH_0, upperBound);
end


