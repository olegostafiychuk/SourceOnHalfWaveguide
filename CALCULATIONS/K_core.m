function K_beta_s_m_talpha_ts_tm =...
          K_core(typeOfCylinder, q, p, psi, m, H2m, dH2m, H1m, dH1m,...
                                tq, tp, tm,...
                 k_0, k, a_0, EE1, GG1, HH1, MU1, EE2, MU2, EE, MU, z0, B_1, B_2, Cm2, Dm2, AE_0, AH_0)

switch(typeOfCylinder)
    
    case 'PerfectConductivity'
        m0 = m;
        q_0 = tq;
        p_0 = tp;
        
        Q = k_0.* a_0 * q;
        H2m  = besselh(m, 2, Q);
        dH2m = (H2m * m)./ Q  - besselh(m + 1, 2, Q);
        H1m  = besselh(m, 1, Q);
        dH1m = (H1m * m)./ Q  - besselh(m + 1, 1, Q);
        
        %%%%%%%% forward wave %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%% first type of waves %%%%%%%%%%%%%%%%%%%%%%%%%%%        
        psi1 = psi;
        
        A = 1./ (k_0 * (q.^2));  
        Ez_q  = Cm2.* (psi1.* H1m + H2m);
        Hz_q  = Dm2.* (-psi1.* H1m + H2m);
        dEz_q = Cm2.* (psi1.* dH1m + dH2m);
        dHz_q = Dm2.* (-psi1.* dH1m + dH2m);
        Ephi_q  = A.* ((-p.* (m./a_0)).* Ez_q + 1i * k_0 * q.* dHz_q);
        Hphi_q  = A.* ((-p.* (m./a_0)).* Hz_q - 1i * k_0 * q.* dEz_q);
    
        Ez_q_0   = AE_0.* besselj(m0,  k_0.* a_0* q_0);
        Hz_q_0   = AH_0.* besselj(m0,  k_0.* a_0* q_0);
        dEz_q_0  = AE_0.* k_0.* (q_0).*((besselj(m0, k_0.* a_0* q_0) * (m0))./ (k_0.* a_0* q_0)  - besselj(m0 + 1,   k_0.* a_0* q_0));
        dHz_q_0  = AH_0.* k_0.* (q_0).*((besselj(m0, k_0.* a_0* q_0) * (m0))./ (k_0.* a_0* q_0)  - besselj(m0 + 1,   k_0.* a_0* q_0));
        A0 = 1./ (k_0 * (1 - p_0.^2));
        Ephi_q_plus0  = A0.* (-p_0.* (m0./a_0).* Ez_q_0 + 1i * dHz_q_0);
        Hphi_q_plus0  = A0.* (-p_0.* (m0./a_0).* Hz_q_0 - 1i * dEz_q_0);  
    
        %%%%%%%%%%%%% exitation coefficient of forward waves %%%%%%%%%%%%%%
        S_plus = (-Ez_q.* Hphi_q_plus0 + Ephi_q.* Hz_q_0 +...
                 Ez_q_0.* Hphi_q - Ephi_q_plus0.* Hz_q).* (a_0./(1i*k_0*(p + p_0)));
             
        K_beta_s_m_talpha_ts_tm = S_plus.* exp(-1i*(p + tp)*k_0*z0);
        
    case 'Isotropic'      
%         [B_1, B_2, Cm2, Dm2]  = coefficientsOfContinuousSpectrum(typeOfCylinder, k_0, k, a_0, EE1, GG1, HH1, MU1, EE, MU, m, p, psi,  1);
        
        % r = a_0;
        %%% структура поля моды
%         H2m  = besselh(m, 2, k_0.* a_0 * q);
%         dH2m = (besselh(m, 2, k_0.* a_0 * q) * m)./ (k_0.* a_0 * q)  - besselh(m + 1, 2, k_0.* a_0 * q);
%         H1m  = besselh(m, 1, k_0.* a_0 * q);
%         dH1m = (besselh(m, 1, k_0.* a_0 * q) * m)./ (k_0.* a_0 * q)  - besselh(m + 1, 1, k_0.* a_0 * q);

        
        A = 1./ (k_0 * (q.^2));  
        Ez_q  = q.* Cm2.* (psi.* H1m + H2m);
        Hz_q  = q.* Dm2.* (-psi.* H1m + H2m);
        Ephi_q  = A.* ((-p.* (m./a_0)).* Ez_q + 1i * k_0 * q.* q.*...
            Dm2.* (-psi.* dH1m + dH2m));
        Hphi_q  = A.* ((-p.* (m./a_0)).* Hz_q - 1i * k_0 * q.* q.*...
            Cm2.* (psi.* dH1m + dH2m));
        
        %%%%% type Of Bessel Beam
        A0 = 1./ (k_0 * (tq.^2));

        Jtmtq    = besselj(tm  , k_0.* a_0* tq);
        Jtmpl1tq = besselj(tm+1, k_0.* a_0* tq);
        Jtm_1tq  = besselj(tm-1, k_0.* a_0* tq);
        Jtmpl2tq = besselj(tm+2, k_0.* a_0* tq);
        
        Ez_tq   = AE_0.* Jtmtq;
        Hz_tq   = AH_0.* Jtmtq;
        dEz_tq  = AE_0.* k_0.* (tq).*((Jtmtq.* (tm))./ (k_0.* a_0* tq)  - Jtmpl1tq);
        dHz_tq  = AH_0.* k_0.* (tq).*((Jtmtq.* (tm))./ (k_0.* a_0* tq)  - Jtmpl1tq);
        
        Ephi_q_plus0  = A0.* ((-tp.* (tm./a_0)).* Ez_tq + 1i * dHz_tq);
        Hphi_q_plus0  = A0.* (-1i * dEz_tq - (tp.* (tm./a_0)).* Hz_tq);
        %     Hrho_q_plus0  = @(r)  A0.* (((tm./r)).* Ez_tq(r) - 1i * tp * dHz_tq(r));
        %     Erho_q_plus0  = @(r)  A0.* (-1i * tp * dEz_tq(r) - (tm./r).* Hz_tq(r));
        
        %%%% из теоремы Гаусса
        S_plus = (-Ez_q.* Hphi_q_plus0 + Ephi_q.* Hz_tq +...
                  Ez_tq.* Hphi_q - Ephi_q_plus0.* Hz_q).* (a_0./(1i*k_0*(p+tp)));
        
        
        
        %     %%% структура поля моды
        q1 = sqrt(EE1*MU1 - p.^2);
        A_1 = 1./ (k_0 * (EE1*MU1 - p.^2));
        
        % %%% структура тестового поля
        q10 = sqrt(EE1*MU1 -  tp.^2);
        A_0 = 1./ (k_0 * (tq.^2));
        
        Jtmq1    = besselj(tm  , k_0.* a_0* q1);
        Jtm_1q1  = besselj(tm-1, k_0.* a_0* q1);
        Jtmpl1q1 = besselj(tm+1, k_0.* a_0* q1);
        Jtmpl2q1 = besselj(tm+2, k_0.* a_0* q1);
        
        I02 = 1/2*((a_0)./ ((k_0*tq).^2 - (k_0*q1).^2)).* (k_0.* tq.* Jtmtq.* Jtm_1q1  - k_0.* q1.* Jtm_1tq.* Jtmq1)-...
            1/2*((a_0)./ ((k_0*tq).^2 - (k_0*q1).^2)).* (k_0.* tq.* Jtmpl2tq.* Jtmpl1q1  - k_0.* q1.* Jtmpl1tq.* Jtmpl2q1);
        I03 = ((a_0)./ ((k_0*tq).^2 - (k_0*q1).^2)).* (k_0.* tq.* Jtmpl2tq.* Jtmpl1q1  - k_0.* q1.* Jtmpl1tq.* Jtmpl2q1);
        I04 =   Jtmtq.* Jtmq1;
        S_minus = (-1)^tm * (A_1.* A_0.* AH_0.* ((tp - p).* B_2).* I02.* (k_0.*tq).* (k_0.*q1)+...
            A_1.* A_0.* AH_0.* ((tp - p).* B_2).* I03.* (k_0.*tq).* (k_0.*q1)+...
            A_1.* A_0.* AH_0.* (1i * (tm).* (EE1 - p.* tp)).* B_1.* I04);
        
        
        S_minus = S_minus -(-1)^tm * A_1.* A_0.* AE_0.* ((((tp*EE1 - p).* B_1)).* (I03.* (k_0.*tq).* (k_0.*q1) + tm * I04)+...
                                            ((1i * (tm).* (p.* tp - 1)).* B_2).* I04);
        S_minus = S_minus.* q1;
                 

        K_beta_s_m_talpha_ts_tm = 2*pi * (S_plus + S_minus).* exp(-1i*(p + tp)*k_0*z0);
        
    case 'Gyrotropic'
        
        [q1, q2, n1, n2, alp1, alp2, bet1, bet2] = term_of_gyrotropic_waveguide(EE1, GG1, HH1, p);
%         [B_1, B_2, Cm2, Dm2]  = coefficientsOfContinuousSpectrum(typeOfCylinder, k_0, k, a_0, EE1, -GG1, HH1, MU1, EE, MU, -m, -p, psi,  1);
        %%%%%% norm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%         % r = a_0;
%         %%% структура поля моды
%         H2m  = besselh(m, 2, k_0.* a_0 * q);
%         dH2m = (besselh(m, 2, k_0.* a_0 * q) * (m))./ (k_0.* a_0 * q)  - besselh(m + 1, 2, k_0.* a_0 * q);
%         H1m  = besselh(m, 1, k_0.* a_0 * q);
%         dH1m = (besselh(m, 1, k_0.* a_0 * q) * (m))./ (k_0.* a_0 * q)  - besselh(m + 1, 1, k_0.* a_0 * q);
        
        A = 1./ (k_0 * (q.^2));  
        Ez_q  = Cm2.* q.* (psi.* H1m + H2m);
        Hz_q  = Dm2.* q.* (-psi.* H1m + H2m);
        Ephi_q  = A.* ((-p.* (m./a_0)).* Ez_q + 1i * k_0 * q.*...
            Dm2.* q.* (-psi.* dH1m + dH2m));
        Hphi_q  = A.* ((-p.* (m./a_0)).* Hz_q - 1i * k_0 * q.*...
            Cm2.* q.* (psi.* dH1m + dH2m));
        
        %%%%% type Of Bessel Beam is H-polarized beam
        A0 = 1./ (k_0 * (tq.^2));
        
        Ez_tq   = AE_0.* besselj(tm,  k_0.* a_0* tq);
        Hz_tq   = AH_0.* besselj(tm,  k_0.* a_0* tq);
        dEz_tq  = AE_0.* k_0.* (tq).*((besselj(tm, k_0.* a_0* tq).* (tm))./ (k_0.* a_0* tq)  - besselj(tm + 1,   k_0.* a_0* tq));
        dHz_tq  = AH_0.* k_0.* (tq).*((besselj(tm, k_0.* a_0* tq).* (tm))./ (k_0.* a_0* tq)  - besselj(tm + 1,   k_0.* a_0* tq));
        
        Ephi_q_plus0  = A0.* ((-tp.* (tm./a_0)).* Ez_tq + 1i * dHz_tq);
        Hphi_q_plus0  = A0.* (-1i * dEz_tq - (tp.* (tm./a_0)).* Hz_tq);
        %     Hrho_q_plus0  = @(r)  A0.* (((tm./r)).* Ez_tq(r) - 1i * tp * dHz_tq(r));
        %     Erho_q_plus0  = @(r)  A0.* (-1i * tp * dEz_tq(r) - (tm./r).* Hz_tq(r));
        
%         %%%% из теоремы Гаусса
%         S_plus = (-1)^m * (-Ez_q.* Hphi_q_plus0 + Ephi_q.* Hz_tq +...
%                            Ez_tq.* Hphi_q - Ephi_q_plus0.* Hz_q).*...
%                           (a_0 * exp(-1i*(-p+tp)*k_0*z)./(1i*k_0*(-p+tp)));
                      
        S_plus = (-Ez_q.* Hphi_q_plus0 + Ephi_q.* Hz_tq +...
                  Ez_tq.* Hphi_q - Ephi_q_plus0.* Hz_q).*...
                  (a_0./(1i*k_0*(p+tp)));    
          
         
%       %     %%% структура поля моды
          JM0    = besselj(m+1, (tq.* a_0).* k_0);                            
          JM1    = besselj(m+1, (q1.* a_0).* k_0);
          JM2    = besselj(m+1, (q2.* a_0).* k_0);
          JMM0   = besselj(m+2, (tq.* a_0).* k_0);
          JMM1   = besselj(m+2, (q1.* a_0).* k_0);
          JMM2   = besselj(m+2, (q2.* a_0).* k_0);
          Jm0    = besselj(m, (tq.* a_0).* k_0);
          Jm1    = besselj(m, (q1.* a_0).* k_0);
          Jm2    = besselj(m, (q2.* a_0).* k_0);
          
             %%%% vary larger JM1 and Jm1
          JMM1(abs(imag(q1.* a_0.* k_0))>300) = besselj(m+2, 1i*300);
          JM1(abs(imag(q1.* a_0.* k_0))>300) = besselj(m+1, 1i*300);
          Jm1(abs(imag(q1.* a_0.* k_0))>300) = besselj(m,   1i*300);
          
          %%%%% added by Oleg Ostafiychuk %%%%%%%%%%%
          JMM2(abs(imag(q2.* a_0.* k_0))>300) = besselj(m+2, 1i*300);
          JM2(abs(imag(q2.* a_0.* k_0))>300) = besselj(m+1, 1i*300);
          Jm2(abs(imag(q2.* a_0.* k_0))>300) = besselj(m,   1i*300);
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
          
          M10 =(a_0./ ((k_0*q1).^2 - (k_0*tq).^2)).*...
              (k_0*q1.* JMM1.* JM0 - k_0*tq.* JMM0.* JM1);
          M20 =(a_0./ ((k_0*q2).^2 - (k_0*tq).^2)).*...
              (k_0*q2.* JMM2.* JM0 - k_0*tq.* JMM0.* JM2);
          
          A_0 = 1./ (k_0 * (tq.^2));
          A_1 = 1./ (k_0 * (GG1^2 - (p.^2 - EE1).^2));
          
           I_into_1 = (-p.* (EE1 - p.^2) - tp* (GG1^2 - EE1.* (EE1 - p.^2)))*...
              (1i/HH1)* AE_0.*...
              ((B_1.* n1.* q1).* (m.* Jm1.* Jm0 + k_0^2 * q1.* tq.* M10)+...
               (B_2.* n2.* q2).* (m.* Jm2.* Jm0 + k_0^2 * q2.* tq.* M20));
          
          I_into_2 = ((EE1 - p.^2).* (p - tp))*...
              AH_0.*...
              ((-B_1.* q1).* (m.* Jm1.* Jm0 + k_0^2 * q1.* tq.* M10)+...
               (-B_2.* q2).* (m.* Jm2.* Jm0 + k_0^2 * q2.* tq.* M20));
           
          I_into_3 = m * (p.* tp - 1).* p.* GG1.*...
                  ((1i./HH1).* (B_1.* n1.* q1.* Jm1 + B_2.* n2.* q2.* Jm2)).*...
                               (AE_0.* Jm0);
           
          I_into_4 = GG1 * (p - tp).*...
              m.* (-B_1.* q1.* Jm1 - B_2.* q2.* Jm2).* (AH_0.* Jm0);
          
          I_into_5 = 1i * (-p.* tp + 1).* (EE1 - p.^2).*...
              m.* (-B_1.* q1.* Jm1 - B_2.* q2.* Jm2).* (AE_0.* Jm0);
          
          I_into_6 = 1i * (-p.* tp.* (EE1 - p.^2) - (GG1^2 - EE1.* (EE1 - p.^2))).*...
              m.* ((1i./HH1).* (B_1.* n1.* q1.* Jm1 + B_2.* n2.* q2.* Jm2)).* (AH_0.* Jm0);

          I_into_7 = 1i * GG1 * (-p.*tp + 1).*...
               AE_0.*...
              ((-B_1.* q1).* (m.* Jm1.* Jm0 + k_0^2 * q1.* tq.* M10)+...
               (-B_2.* q2).* (m.* Jm2.* Jm0 + k_0^2 * q2.* tq.* M20));
           
           I_into_8 = 1i * GG1 * p.* (-tp + p).*...
              (1i/HH1)* AH_0.*...
              ((B_1.* n1.* q1).* (m.* Jm1.* Jm0 + k_0^2 * q1.* tq.* M10)+...
               (B_2.* n2.* q2).* (m.* Jm2.* Jm0 + k_0^2 * q2.* tq.* M20));
          
          S_minus = (-1)^m * (A_1.* A_0).* (I_into_1 + I_into_2 + I_into_3 +...
                                            I_into_4 + I_into_5 + I_into_6 +...
                                            I_into_7 + I_into_8);

         K_beta_s_m_talpha_ts_tm = 2*pi * (S_plus + S_minus).* exp(-1i*(p + tp)*k_0*z0);
     
end


