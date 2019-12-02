function N_n = Norm_of_descreteMode(typeOfCylinder, q, p, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m)

c = 3e10;   % velocity of light

  switch(typeOfCylinder)
      case 'Isotropic'
          m = abs(m);
          %%% analytical formula
          [B_1_mn,B_2_mn,Cm2_mn, Dm2_mn] = coefficientsOfContinuousSpectrum(typeOfCylinder, k_0, k_0, a_0, EE1,0,0, MU1, EE, MU, -m, (-p), q, 0, 1);
          [B_1_pl,B_2_pl,Cm2_pl, Dm2_pl] = coefficientsOfContinuousSpectrum(typeOfCylinder, k_0, k_0, a_0, EE1,0,0, MU1, EE, MU, m, (p), q, 0, 1);
          %%% norm
          q1 = sqrt(EE1*MU1/(EE*MU) -  p.^2);
          A_1 = 1./ (k_0 * (EE1*MU1 - EE*MU * p.^2));       
          A = 1./ (k_0 * (1 - p.^2)); 
          
          G_Jm_a0 = - 1/2 * (besselj(0, k_0*q1*a_0)).^2;
          G_Jm_0  = - 1/2 * (besselj(0, k_0*q1*0)).^2;    
          G_Hm2_a0 = - 1/2 * (besselh(0,2, k_0*q*a_0)).^2;
          for n=0:(m-1)
              G_Jm_a0  =  G_Jm_a0  - (besselj(m-n, k_0*q1*a_0)).^2;
              G_Hm2_a0 =  G_Hm2_a0 - (besselh(m-n, 2, k_0*q*a_0)).^2;
          end
          
          G_Jm_1__a0 = - 1/2 * (besselj(0, k_0*q1*a_0)).^2;
          G_Jm_1__0  = - 1/2 * (besselj(0, k_0*q1*0)).^2;
          G_Hm2_1__a0 = - 1/2 * (besselh(0,2, k_0*q*a_0)).^2;    
          for n=0:(m-2)
              G_Jm_1__a0  =  G_Jm_1__a0  - (besselj(m-1-n, k_0*q1*a_0)).^2;
              G_Hm2_1__a0 =  G_Hm2_1__a0 - (besselh(m-1-n, 2, k_0*q*a_0)).^2;
          end
          
          F_Jm_a0 = ((k_0*q1*a_0).^2 / 2).* ((besselj(m+1, k_0*q1*a_0)).^2  - besselj(m,   k_0*q1*a_0).* besselj(m+2,  k_0*q1*a_0));
          F_Hm2_a0 =((k_0*q *a_0).^2 / 2).* ((besselh(m+1,2, k_0*q*a_0)).^2 - besselh(m,2, k_0*q *a_0).* besselh(m+2,2, k_0*q *a_0));
          
          I_into = (-1)^m * (A_1.^2).* (-2 * p)* (EE1 * B_1_mn.* B_1_pl - B_2_mn.* B_2_pl).*...
              (-m * (G_Jm_a0 - G_Jm_0) + m * (G_Jm_1__a0 - G_Jm_1__0) + F_Jm_a0);
          I_into = I_into + (-1)^m * (A_1.^2) * m * 1i * (p.^2 + EE1*MU1) *...
              ((B_1_mn.* B_2_pl + B_2_mn.* B_1_pl).* (besselj(m, k_0*q1*a_0)).^2);%%% добавление полного дифференциала
          
          
          I_out  =  (-1)^m * (A.^2).* (-2 * p)*(Cm2_mn.* Cm2_pl - Dm2_mn.* Dm2_pl).*...
              (-m * G_Hm2_a0 + m * G_Hm2_1__a0 + F_Hm2_a0);
          I_out  = I_out + (-1)^m * (A.^2) * m * 1i * (p.^2 + 1) *...
              ((Cm2_mn.* Dm2_pl + Cm2_pl.* Dm2_mn).* (besselh(m,2, k_0*q*a_0)).^2);%%% добавление полного дифференциала
          
          I_into= q1.^2.* I_into; %%% домножаем на q1^2 из-за представления
          I_out = q.^2.*  I_out; %%% домножаем на q^2 из-за представления
          
          N_n = (I_into - I_out) * c / 2;
          
          
      case 'Gyrotropic'          
             [q1, q2, n1, n2, alp1, alp2, bet1, bet2] = term_of_gyrotropic_waveguide(EE1, GG1, HH1, p);
             [B_1,B_2,Cm2, Dm2] = coefficientsOfContinuousSpectrum('Gyrotropic', k_0, k_0, a_0, EE1, GG1,HH1, MU1, EE, MU,  m, (p), q, 0, 1);
             
             JM1    = besselj(m+1, (q1.* a_0).* k_0);
             JM2    = besselj(m+1, (q2.* a_0).* k_0);
             JMM1    = besselj(m+2, (q1.* a_0).* k_0);
             JMM2    = besselj(m+2, (q2.* a_0).* k_0);
             Jm1    = besselj(m, (q1.* a_0).* k_0);
             Jm2    = besselj(m, (q2.* a_0).* k_0);
             Jm1_Q1 = Jm1./ ((q1.* a_0).* k_0);
             Jm2_Q2 = Jm2./ ((q2.* a_0).* k_0);
             
             F1 = ((k_0*q1*a_0).^2 / 2).* (JM1.^2  - Jm1.*  JMM1);
             F2 = ((k_0*q2*a_0).^2 / 2).* (JM2.^2  - Jm2.*  JMM2);
             M12 =(a_0./ ((k_0*q1).^2 - (k_0*q2).^2)).*...
                 (k_0*q1.* JMM1.* JM2 - k_0*q2.* JMM2.* JM1);
    
             
             A_1 = 1./ (k_0 * (GG1^2 - (p.^2 - EE1).^2));
             
             I_into_1 = (-p.^3 * GG1^2 + p.* (EE1 - p.^2).* (GG1^2 - EE1.* (EE1 - p.^2)))*...
                 (-1/HH1^2)*...
                 ((B_1.* n1.* q1).^2 * (m.* Jm1.^2 + F1)+...
                 (B_2.* n2.* q2).^2 * (m.* Jm2.^2 + F2)+...
                 2*(B_1.* n1.* q1).* (B_2.* n2.* q2).* (m.* Jm1.* Jm2 + k_0^2 * q1.* q2.* M12));
             
             I_into_2 = (p * GG1^2 + p.* (EE1 - p.^2).^2)*...
                 ((B_1.* q1).^2 * (m.* Jm1.^2 + F1)+...
                 (B_2.* q2).^2 * (m.* Jm2.^2 + F2)+...
                 2*(B_1.* q1).* (B_2.* q2).* (m.* Jm1.* Jm2 + k_0^2 * q1.* q2.* M12));
             
             I_into_3 = 1i* (-(EE1 - p.^2).* (GG1^2 - EE1.* (EE1 - p.^2)) +...
                 p.^2 * (EE1 - p.^2).^2 +...
                 2 * p.^2* GG1^2).*...
                 m.* (1i./HH1).* (B_1.* n1.* q1.* Jm1 +  B_2.* n2.* q2.* Jm2).*...
                 (-B_1.* q1.* Jm1 - B_2.* q2.* Jm2);
             
             I_into_4 = (p.*GG1.*(GG1^2 - EE1.* (EE1 - p.^2)) - p.^3 * GG1 * (EE1 - p.^2))*...
                 m * (-1./HH1^2) * ((B_1.* n1.* q1.* Jm1 + B_2.* n2.* q2.* Jm2)).^2;
             
             I_into_5 = 2*(p.*GG1.* (EE1 - p.^2))*...
                 m * (B_1.* q1.* Jm1 + B_2.* q2.* Jm2).^2;
             
             I_into_6 = -(-1i.*GG1.*(GG1^2 - EE1.* (EE1 - p.^2))+3*1i* p.^2 * GG1 * (EE1 - p.^2)).*...
                 (1i./HH1).* ((B_1.* q1).^2.*n1.* (m.* Jm1.^2 + F1)+...
                 (B_2.* q2).^2.*n2.* (m.* Jm2.^2 + F2)+...
                 (B_1.* q1).* (B_2.* q2).* (n1 + n2).* (m.* Jm1.* Jm2 + k_0^2 * q1.* q2.* M12));
             
             I_into = (-1)^m * (A_1.^2).* 2.* (I_into_1 + I_into_2 + I_into_3 +...
                 I_into_4 + I_into_5 + I_into_6);
             

          m_Old = m;
          m = abs(m);
          A = 1./ (k_0 * (1 - p.^2));
          G_Hm2_a0 = - 1/2 * (besselh(0,2, k_0*q*a_0)).^2;
          for n=0:(m-1)
              G_Hm2_a0 =  G_Hm2_a0 - (besselh(m-n, 2, k_0*q*a_0)).^2;
          end
          
          G_Hm2_1__a0 = - 1/2 * (besselh(0,2, k_0*q*a_0)).^2;
          for n=0:(m-2)
              G_Hm2_1__a0 =  G_Hm2_1__a0 - (besselh(m-1-n, 2, k_0*q*a_0)).^2;
          end
          
          F_Hm2_a0 =((k_0*q *a_0).^2 / 2).* ((besselh(m+1,2, k_0*q*a_0)).^2 - besselh(m,2, k_0*q *a_0).* besselh(m+2,2, k_0*q *a_0));
          
          I_out  =  (-1)^m * (A.^2).* (-2 * p)*(Cm2.* Cm2 - Dm2.* Dm2).*...
              (-m * G_Hm2_a0 + m * G_Hm2_1__a0 + F_Hm2_a0);
          I_out  = I_out + (-1)^m * (A.^2) * m * 1i * (p.^2 + 1) *...
              ((Cm2.* Dm2 + Cm2.* Dm2).* (besselh(m,2, k_0*q*a_0)).^2)*...
              (sign(m_Old));
          
          I_out = q^2 * I_out; %%% домножаем на q^2 из-за представления
          
          N_n = (I_into - I_out) * c / 2;
  end





