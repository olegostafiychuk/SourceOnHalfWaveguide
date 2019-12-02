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
          coeffsOfFieldAndExcitationCoeffs_of_Continues_discreteRepr(typeOfCylinder, q, q_0, k_0, k, a_0, EE1, GG1, HH1, MU1, EE2, MU2, EE, MU, m, z, AE_0, AH_0)

  p = sqrt(1-q.^2);
  p = real(p) - 1i * abs(imag(p));
  
  p_0 = sqrt(1-q_0.^2);
  p_0 = real(p_0) - 1i * abs(imag(p_0));
  
  m0 = m;
  
switch(typeOfCylinder)
    case 'PerfectConductivity'
        
        m = -m;
        p = -p;
        
        Q = k_0.* a_0 * q;
        H2m  = besselh(m, 2, Q);
        dH2m = (H2m * m)./ Q  - besselh(m + 1, 2, Q);
        H1m  = besselh(m, 1, Q);
        dH1m = (H1m * m)./ Q  - besselh(m + 1, 1, Q);
        
        %%%%%%%% forward wave %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%% first type of waves %%%%%%%%%%%%%%%%%%%%%%%%%%%        
        psi1 = -besselh(m,2, k_0 * q * a_0)./ besselh(m,1, k_0 * q * a_0);
        Dm2_forward_1 = 0;
        Cm2_forward_1 = 1;
        N_p_1 = (-1)^m * 8 * (p./ q.^3 / k_0).* (Cm2_forward_1.^2 + Dm2_forward_1.^2).* psi1;
        
        A = 1./ (k_0 * (1 - p.^2));  
        Ez_q  = Cm2_forward_1.* (psi1.* H1m + H2m);
        Hz_q  = Dm2_forward_1.* (-psi1.* H1m + H2m);
        dEz_q = Cm2_forward_1.* (psi1.* dH1m + dH2m);
        dHz_q = Dm2_forward_1.* (-psi1.* dH1m + dH2m);
        Ephi_q  = A.* ((-p.* (m./a_0)).* Ez_q + 1i * k_0 * q.* dHz_q);
        Hphi_q  = A.* ((-p.* (m./a_0)).* Hz_q - 1i * k_0 * q.* dEz_q);
    
        Ez_q_0   = AE_0.* besselj(m0,  k_0.* a_0* q_0);
        Hz_q_0   = AH_0.* besselj(m0,  k_0.* a_0* q_0);
        dEz_q_0  = AE_0.* k_0.* (q_0) *((besselj(m0, k_0.* a_0* q_0) * (m0))./ (k_0.* a_0* q_0)  - besselj(m0 + 1,   k_0.* a_0* q_0));
        dHz_q_0  = AH_0.* k_0.* (q_0) *((besselj(m0, k_0.* a_0* q_0) * (m0))./ (k_0.* a_0* q_0)  - besselj(m0 + 1,   k_0.* a_0* q_0));
        A0 = 1./ (k_0 * (1 - p_0.^2));
        Ephi_q_plus0  = A0.* (-p_0.* (m0./a_0).* Ez_q_0 + 1i * dHz_q_0);
        Hphi_q_plus0  = A0.* (-p_0.* (m0./a_0).* Hz_q_0 - 1i * dEz_q_0);
%         Hrho_q_plus0  = @(r)  A0.* (((m0./a_0)).* Ez_q_0 - 1i * p_0 * dHz_q_0);
%         Erho_q_plus0  = @(r)  A0.* (-1i * p_0 * dEz_q_0 - (m0./a_0).* Hz_q_0);    
    
        %%%%%%%%%%%%% exitation coefficient of forward waves %%%%%%%%%%%%%%
        S_plus = (-Ez_q.* Hphi_q_plus0 + Ephi_q.* Hz_q_0 +...
                 Ez_q_0.* Hphi_q - Ephi_q_plus0.* Hz_q).* (a_0 * exp(-1i*(p+p_0)*k_0*z)./(1i*k_0*(p + p_0)));
        a_p_field_1_forw = S_plus./ N_p_1;
        
        
        
        %%%%%%%%%%%%% exitation coefficient of backward waves %%%%%%%%%%%%%%
        Ephi_q  = A.* (( p.* (m./a_0)).* Ez_q + 1i * k_0 * q.* dHz_q);
        Hphi_q  = A.* (( p.* (m./a_0)).* Hz_q - 1i * k_0 * q.* dEz_q);
        S_plus = (-Ez_q.* Hphi_q_plus0 + Ephi_q.* Hz_q_0 +...
                 Ez_q_0.* Hphi_q - Ephi_q_plus0.* Hz_q).* (a_0 * exp(-1i*(-p+p_0)*k_0*z)./(1i*k_0*(-p + p_0)));
        a_p_field_1_back = -S_plus./ N_p_1;
        
        
        
        %%%%%%%%%%%%%%%%%% second type of waves %%%%%%%%%%%%%%%%%%%%%%%%%%%
        psi2 =  (besselh(m,2, k_0 * q * a_0) * m./ (k_0 * q * a_0) - besselh(m+1,2, k_0 * q * a_0))./...
                (besselh(m,1, k_0 * q * a_0) * m./ (k_0 * q * a_0) - besselh(m+1,1, k_0 * q * a_0));
        Dm2_forward_2 = 1;
        Cm2_forward_2 = 0;
        N_p_2 = (-1)^m * 8 * (p./ q.^3 / k_0).* (Cm2_forward_2.^2 + Dm2_forward_2.^2).* psi2;
        
        A = 1./ (k_0 * (1 - p.^2));  
        Ez_q  = Cm2_forward_2.* (psi2.* H1m + H2m);
        Hz_q  = Dm2_forward_2.* (-psi2.* H1m + H2m);
        dEz_q = Cm2_forward_2.* (psi2.* dH1m + dH2m);
        dHz_q = Dm2_forward_2.* (-psi2.* dH1m + dH2m);
        Ephi_q  = A.* ((-p.* (m./a_0)).* Ez_q + 1i * k_0 * q.* dHz_q);
        Hphi_q  = A.* ((-p.* (m./a_0)).* Hz_q - 1i * k_0 * q.* dEz_q);
        
        %%%%%%%%%%%%% exitation coefficient of forward waves %%%%%%%%%%%%%%
        S_plus = (-Ez_q.* Hphi_q_plus0 + Ephi_q.* Hz_q_0 +...
                 Ez_q_0.* Hphi_q - Ephi_q_plus0.* Hz_q).* (a_0 * exp(-1i*(p+p_0)*k_0*z)./(1i*k_0*(p + p_0)));
        a_p_field_2_forw = S_plus./ N_p_2;
        
        %%%%%%%%%%%%% exitation coefficient of backward waves %%%%%%%%%%%%%%
        Ephi_q  = A.* (( p.* (m./a_0)).* Ez_q + 1i * k_0 * q.* dHz_q);
        Hphi_q  = A.* (( p.* (m./a_0)).* Hz_q - 1i * k_0 * q.* dEz_q);
        S_plus = (-Ez_q.* Hphi_q_plus0 + Ephi_q.* Hz_q_0 +...
                 Ez_q_0.* Hphi_q - Ephi_q_plus0.* Hz_q).* (a_0 * exp(-1i*(-p+p_0)*k_0*z)./(1i*k_0*(-p + p_0)));
        a_p_field_2_back = -S_plus./ N_p_2;
        
        psi_forward_1  = psi1;
        psi_backward_1 = psi1;
        psi_forward_2  = psi2;
        psi_backward_2 = psi2;
        psi_forward_1_transp  = [];
        psi_backward_1_transp = [];
        psi_forward_2_transp  = [];
        psi_backward_2_transp = [];
        
        B_1_forward_1 = 0;B_1_forward_2 = 0;
        B_2_forward_1 = 0;B_2_forward_2 = 0;
        B_1_backward_1 = 0;B_1_backward_2 = 0;
        B_2_backward_1 = 0;B_2_backward_2 = 0;
        Cm2_backward_1 = Cm2_forward_1;
        Cm2_backward_2 = Cm2_forward_2;
        Dm2_backward_1 = Dm2_forward_1;
        Dm2_backward_2 = Dm2_forward_2;
        
        
    case 'Isotropic'      
        %%% вычисляем первый тип волн непрерывного спектра
       
        %%% first type of waves
        psi_forward_1  = psi1_q(k_0, k, a_0, EE1, MU1, EE, MU, m,   p, q);
        psi_backward_1  = psi1_q(k_0, k, a_0, EE1, MU1, EE, MU, m,  -p, q);
        [B_1_forward_1,B_2_forward_1,Cm2_forward_1,  Dm2_forward_1]  = coefficientsOfContinuousSpectrum(typeOfCylinder,k_0, k, a_0, EE1,0,0, MU1, EE, MU, m,   p, q,   psi_forward_1,  1);
        [B_1_backward_1,B_2_backward_1,Cm2_backward_1, Dm2_backward_1] = coefficientsOfContinuousSpectrum(typeOfCylinder,k_0, k, a_0, EE1,0,0, MU1, EE, MU, m, (-p), q, psi_backward_1, 1);
%         [B_1_back_PandM_1,B_2_back_PandM_1,Cm2_back_PandM_1, Dm2_back_PandM_1] = coefficientsOfContinuousSpectrum(typeOfCylinder,k_0, k, a_0, EE1, MU1, EE, MU, -m, (-p), q, psi_forward_1, 1);
        
        %%% second type of waves
        psi_forward_2  = psi2_q(k_0, k, a_0, EE1, MU1, EE, MU, m,   p, q);
        psi_backward_2  = psi2_q(k_0, k, a_0, EE1, MU1, EE, MU, m,  -p, q);
        [B_1_forward_2,B_2_forward_2,Cm2_forward_2,  Dm2_forward_2]  = coefficientsOfContinuousSpectrum_forPsi2(typeOfCylinder,k_0, k, a_0, EE1,0,0, MU1, EE, MU, m,   p, q,   psi_forward_2,  1);
        [B_1_backward_2,B_2_backward_2,Cm2_backward_2, Dm2_backward_2] = coefficientsOfContinuousSpectrum_forPsi2(typeOfCylinder,k_0, k, a_0, EE1,0,0, MU1, EE, MU, m, (-p), q, psi_backward_2, 1);
%         [B_1_back_PandM_2,B_2_back_PandM_2,Cm2_back_PandM_2, Dm2_back_PandM_2] = coefficientsOfContinuousSpectrum_forPsi2(typeOfCylinder,k_0, k, a_0, EE1, MU1, EE, MU, -m, (-p), q, psi_forward_2, 1);
        
        psi_forward_1_transp  = [];
        psi_backward_1_transp = [];
        psi_forward_2_transp  = [];
        psi_backward_2_transp = [];

        %%%%% forward waves first type of waves %%%%%%%%%%%%%%%%%%%%%%%%%%%
        a_p_field_1_forw = a_p_field__ForBesselBeam(typeOfCylinder, q, q_0,  p, k_0, k, a_0, EE1, 0, 0, MU1, EE2, MU2, EE, MU, m, z,psi_forward_1, AE_0, AH_0);
        %%%%% backward waves first type of waves %%%%%%%%%%%%%%%%%%%%%%%%%%%
        a_p_field_1_back = a_p_field__ForBesselBeam(typeOfCylinder, q, q_0, -p, k_0, k, a_0, EE1, 0, 0, MU1, EE2, MU2, EE, MU, m, z,psi_backward_1, AE_0, AH_0);
                                          
        %%%%% forward waves second type of waves %%%%%%%%%%%%%%%%%%%%%%%%%%%
        a_p_field_2_forw = a_p_field__ForBesselBeam(typeOfCylinder, q, q_0,  p, k_0, k, a_0, EE1, 0, 0, MU1, EE2, MU2, EE, MU, m, z,psi_forward_2, AE_0, AH_0);
        %%%%% backward waves second type of waves %%%%%%%%%%%%%%%%%%%%%%%%%%%
        a_p_field_2_back = a_p_field__ForBesselBeam(typeOfCylinder, q, q_0, -p, k_0, k, a_0, EE1, 0, 0, MU1, EE2, MU2, EE, MU, m, z,psi_backward_2, AE_0, AH_0);
      
        
        
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
        [B_1_forward_2,B_2_forward_2,Cm2_forward_2,  Dm2_forward_2]  = coefficientsOfContinuousSpectrum_forPsi2(typeOfCylinder,k_0, k, a_0, EE1, GG1, HH1, MU1, EE, MU, m,   p, q,  psi_forward_2,  1);
        [B_1_backward_2,B_2_backward_2,Cm2_backward_2, Dm2_backward_2] = coefficientsOfContinuousSpectrum_forPsi2(typeOfCylinder,k_0, k, a_0, EE1, GG1, HH1, MU1, EE, MU, m, (-p), q, psi_backward_2, 1);
%         [B_1_back_PandM_2,B_2_back_PandM_2,Cm2_back_PandM_2, Dm2_back_PandM_2] = coefficientsOfContinuousSpectrum_forPsi2(typeOfCylinder,k_0, k, a_0, EE1, MU1, EE, MU, -m, (-p), q, psi_forward_2, 1);


        %%%%% forward waves first type of waves %%%%%%%%%%%%%%%%%%%%%%%%%%%
        a_p_field_1_forw = a_p_field__ForBesselBeam(typeOfCylinder, q, q_0,  p, k_0, k, a_0, EE1, GG1, HH1,...
            MU1, EE2, MU2, EE, MU, m, z,psi_backward_1_transp, AE_0, AH_0);
        %%%%% backward waves first type of waves %%%%%%%%%%%%%%%%%%%%%%%%%%%
        a_p_field_1_back = a_p_field__ForBesselBeam(typeOfCylinder, q, q_0, -p, k_0, k, a_0, EE1, GG1, HH1,...
            MU1, EE2, MU2, EE, MU, m, z,psi_forward_1_transp, AE_0, AH_0);
                                          
        %%%%% forward waves second type of waves %%%%%%%%%%%%%%%%%%%%%%%%%%%
        a_p_field_2_forw = a_p_field__ForBesselBeam(typeOfCylinder, q, q_0,  p, k_0, k, a_0, EE1, GG1, HH1,...
            MU1, EE2, MU2, EE, MU, m, z, psi_backward_2_transp, AE_0, AH_0);
        %%%%% backward waves second type of waves %%%%%%%%%%%%%%%%%%%%%%%%%%%
        a_p_field_2_back = a_p_field__ForBesselBeam(typeOfCylinder, q, q_0, -p, k_0, k, a_0, EE1, GG1, HH1,...
            MU1, EE2, MU2, EE, MU, m, z, psi_forward_2_transp, AE_0, AH_0);
end


