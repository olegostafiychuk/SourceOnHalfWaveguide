%%%%% расчёт коэффициента для волн непрерывного спектра в праямом
%%%%% направлении (в сторону среды 2)
%%%%% второе приближение в нахождении коэффициента aplus
function aplus2 = a_p_field_descreteMode__ForBesselBeam(typeOfCylinder, q, q_0, p, k_0, k,   a_0, EE1, GG1, HH1, MU1, EE2, MU2, EE, MU, m, m0, z, AE_0, AH_0)

p_0 = sqrt(1-q_0.^2);
p_0 = real(p_0) - 1i * abs(imag(p_0));

c = 3e10;   % velocity of light

m = -m;
p = -p;
GG1 = -GG1;

  switch(typeOfCylinder)
      case 'Isotropic'  
          [B_1, B_2, Cm2, Dm2]  = coefficientsOfContinuousSpectrum(typeOfCylinder, k_0, k_0, a_0, EE1,0,0, MU1, EE, MU,  m,   p, q,  0, 1);
          
          
          %%%%%%%%%%%%%%%%% test %%%%%%%%%%%%%%%%%%%
          % [B_1,B_2,Cm2, Dm2] = coefficientsOfContinuousSpectrum(k_0, k_0, a_0, EE1, MU1, EE, MU, m, (p), 0, 1);
          % [B_01, B_02, DE, DH] =coefficientsOfContinuousSpectrum(k_0, k_0, a_0, EE1, MU1, EE, MU, m0,   p_0,  0, 1);%%% коэффициенты для тестового поля
          
          %%%% inner integral
          q1 = sqrt(EE1*MU1 -  p.^2);
          q10 = sqrt(EE1*MU1 -  p_0.^2);
          A_1 = 1./ (k_0 * (EE1*MU1 - EE*MU * p.^2));
          %     A_10 = 1./ (k_0 * (EE1*MU1 - EE*MU * p_0.^2)); 
          %     A = 1./ (k_0 * (1 - p.^2)); 
          
          Ez_q   = @(r) q1.* B_1.* besselj(m, k_0.* r* q1);
          Hz_q   = @(r) q1.* B_2.* besselj(m, k_0.* r* q1);
          dEz_q  = @(r) q1.* B_1.* k_0.* q1 *((besselj(m, k_0.* r* q1) * m)./ (k_0.* r* q1)  - besselj(m + 1, k_0.* r* q1));
          dHz_q  = @(r) q1.* B_2.* k_0.* q1 *((besselj(m, k_0.* r* q1) * m)./ (k_0.* r* q1)  - besselj(m + 1, k_0.* r* q1));
          
          Ephi_q_minus  = @(r)  A_1.* (( -p.* (m./r)).* Ez_q(r) + 1i * dHz_q(r));
          Hrho_q_minus  = @(r)  A_1.* ((EE1 * (m./r)).* Ez_q(r) - 1i * p * dHz_q(r));
          Hphi_q_minus  = @(r)  A_1.* (-1i * EE1 * dEz_q(r) - (p.* (m./r)).* Hz_q(r));
          Erho_q_minus  = @(r)  A_1.* (-1i * p * dEz_q(r) - (m./r).* Hz_q(r));
          
          %%%%% type Of Bessel Beam is H-polarized beam
          %   Hz_inc    = AE_0.* (besselj(m0, k_0.* a_0 * q_0));                       
          %   Ez_inc    = 0;
          %   Hphi_inc  = AE_0.* ((k_0 * q_0.^2)^(-1)).* (- (p_0.* (m0./a_0)).*  (besselj(m0, k_0.* a_0 * q_0)));
          %   Ephi_inc  = AE_0.* ((k_0 * q_0.^2)^(-1)).* (1i * MU * k_0 * (q_0).*...
          %       (besselj(m0-1, k_0.* a_0 * q_0) - besselj(m0+1, k_0.* a_0 * q_0))/2);
          %   
          %   [B_01,B_02,DE,DH] = singleCylinderCoefficients_BesselBeam_forTest(k_0 * c, a_0, EE1, MU1, EE, MU,  c, 1, m0, p_0, q_0, 0, Ez_inc, Hz_inc, Ephi_inc, Hphi_inc);
          % %   AE_0 = 0;
          A = 1./ (k_0 * (1 - p.^2));
          A0 = 1./ (k_0 * (1 - p_0.^2));
          
          Ez_q0   = @(r) AE_0.* besselj(m0,  k_0.* r* q_0);
          Hz_q0   = @(r) AH_0.* besselj(m0,  k_0.* r* q_0);
          dEz_q0  = @(r) AE_0.* k_0.* (q_0) *((besselj(m0,   k_0.* r* q_0) * (m0))./ (k_0.* r* q_0)  - besselj(m0 + 1,   k_0.* r* q_0));
          dHz_q0  = @(r) AH_0.* k_0.* (q_0) *((besselj(m0,   k_0.* r* q_0) * (m0))./ (k_0.* r* q_0)  - besselj(m0 + 1,   k_0.* r* q_0));
          
          Ephi_q_plus0  = @(r)  A0.* ((-p_0.* (m0./r)).* Ez_q0(r) + 1i * dHz_q0(r));
          Hphi_q_plus0  = @(r)  A0.* (-1i * dEz_q0(r) - (p_0.* (m0./r)).* Hz_q0(r));
          Hrho_q_plus0  = @(r)  A0.* (((m0./r)).* Ez_q0(r) - 1i * p_0 * dHz_q0(r));
          Erho_q_plus0  = @(r)  A0.* (-1i * p_0 * dEz_q0(r) - (m0./r).* Hz_q0(r));
          
          
          S_minus = integral(@(r) (Erho_q_plus0(r).* Hphi_q_minus(r) - Ephi_q_plus0(r).* Hrho_q_minus(r)-...
              Erho_q_minus(r).* Hphi_q_plus0(r) + Ephi_q_minus(r).* Hrho_q_plus0(r)).* r, 0, a_0);
          % %     
          % %     
          % %     Jm0  = @(r) besselj(m0,  k_0.* r* q_0);
          % %     dJm0 = @(r) k_0.* q_0 *((besselj(m0,   k_0.* r* q_0) * (m0))./ (k_0.* r* q_0)  - besselj(m0 + 1,   k_0.* r* q_0));
          % %     Jm   = @(r) besselj(m0,  k_0.* r* q1);
          % %     dJm  = @(r) k_0.* q1 *((besselj(m0,   k_0.* r* q1) * (m0))./ (k_0.* r* q1)  - besselj(m0 + 1,   k_0.* r* q1));
          %     
          % %     S_minus = (-1)^m0 * integral(@(r) A_1.* A0.* AE_0 * ((1i * (m0./r).* (EE1 - p.* p_0)).* B_1.* (Jm0(r).* dJm(r) + dJm0(r).* Jm(r)) +...
          % %         ((p_0 - p) * B_2).* ((m0^2./ r.^2).* Jm0(r).* Jm(r) + dJm0(r).* dJm(r))).* r, 0, a_0)
          % %     S_minus = (-1)^m0 * (integral(@(r) A_1.* A0.* AE_0 * (((p_0 - p) * B_2).* ((m0^2./ r.^2).* Jm0(r).* Jm(r) + dJm0(r).* dJm(r))).* r, 0, a_0)+...
          % %         A_1.* A0.* AE_0 * (1i * (m0).* (EE1 - p.* p_0)).* B_1.* (Jm0(a_0).* Jm(a_0)))
          % %     S_minus = (-1)^m0 * (integral(@(r) A_1.* A0.* AE_0 * (((p_0 - p) * B_2).*...
          % %         ((m0./ r).* (besselj(m0,   k_0.* r* q_0).* besselj(m0-1, k_0.* r* q1)*k_0.*q1 -...
          % %                      besselj(m0+1, k_0.* r* q_0).* besselj(m0,   k_0.* r* q1)*k_0.*q_0) +...
          % %                      besselj(m0+1, k_0.* r* q_0).* besselj(m0+1, k_0.* r* q1) * (k_0.*q_0) * (k_0.*q1))).* r, 0, a_0)+...
          % %         A_1.* A0.* AE_0 * (1i * (m0).* (EE1 - p.* p_0)).* B_1.* (Jm0(a_0).* Jm(a_0)))
          % 
          % I02 = 1/2*((a_0)./ ((k_0*q_0).^2 - (k_0*q1).^2)).* (k_0.* q_0.* besselj(m0, k_0.* a_0* q_0).* besselj(m0-1, k_0.* a_0* q1)  - k_0.* q1.* besselj(m0-1, k_0.* a_0* q_0).* besselj(m0, k_0.* a_0* q1))-...
          %       1/2*((a_0)./ ((k_0*q_0).^2 - (k_0*q1).^2)).* (k_0.* q_0.* besselj(m0+2, k_0.* a_0* q_0).* besselj(m0+1, k_0.* a_0* q1)  - k_0.* q1.* besselj(m0+1, k_0.* a_0* q_0).* besselj(m0+2, k_0.* a_0* q1));
          % I03 = ((a_0)./ ((k_0*q_0).^2 - (k_0*q1).^2)).* (k_0.* q_0.* besselj(m0+2, k_0.* a_0* q_0).* besselj(m0+1, k_0.* a_0* q1)  - k_0.* q1.* besselj(m0+1, k_0.* a_0* q_0).* besselj(m0+2, k_0.* a_0* q1));
          % I04 =   besselj(m0, k_0.* a_0* q_0).* besselj(m0, k_0.* a_0* q1);
          % S_minus = (-1)^m0 * (A_1.* A0.* AH_0 * ((p_0 - p) * B_2).* I02 * (k_0.*q_0) * (k_0.*q1)+...
          %                      A_1.* A0.* AH_0 * ((p_0 - p) * B_2) * I03 * (k_0.*q_0) * (k_0.*q1)+...
          %                      A_1.* A0.* AH_0 * (1i * (m0).* (EE1 - p.* p_0)).* B_1.* I04)
          
          
          %     S_minus = (-Ez_q(a_0).* Hphi_q_plus0(a_0) + Ephi_q_minus(a_0).* Hz_q0(a_0) +...
          %                Ez_q0(a_0).* Hphi_q_minus(a_0) - Ephi_q_plus0(a_0).* Hz_q(a_0)).* (a_0./(1i*k_0*(p + p_0)));
          
          %%%% outer integral   
          Ez_q   = @(r) q.* Cm2.* besselh(m,2, k_0.* r* q);
          Hz_q   = @(r) q.* Dm2.* besselh(m,2, k_0.* r* q);
          dEz_q  = @(r) q.* Cm2.* k_0.* q *((besselh(m,2, k_0.* r* q) * m)./ (k_0.* r* q)  - besselh(m + 1,2, k_0.* r* q));
          dHz_q  = @(r) q.* Dm2.* k_0.* q *((besselh(m,2, k_0.* r* q) * m)./ (k_0.* r* q)  - besselh(m + 1,2, k_0.* r* q));
          
          Ephi_q_out_minus  = @(r)  A.* ((-p.* (m./r)).* Ez_q(r) + 1i * dHz_q(r));
          %     Hrho_q_out_minus  = @(r)  A.* (((m./r)).* Ez_q(r) - 1i * p * dHz_q(r));
          Hphi_q_out_minus  = @(r)  A.* (-1i * dEz_q(r) - (p.* (m./r)).* Hz_q(r));
          %     Erho_q_out_minus  = @(r)  A.* (-1i * p * dEz_q(r) - (m./r).* Hz_q(r));
          
          
          
          %%%%% the Bessel Beam
          Ez_q0   = @(r) AE_0.* besselj(m0,  k_0.* r* q_0);
          Hz_q0   = @(r) AH_0.* besselj(m0,  k_0.* r* q_0);
          dEz_q0  = @(r) AE_0.* k_0.* (q_0) *((besselj(m0,   k_0.* r* q_0) * (m0))./ (k_0.* r* q_0)  - besselj(m0 + 1,   k_0.* r* q_0));
          dHz_q0  = @(r) AH_0.* k_0.* (q_0) *((besselj(m0,   k_0.* r* q_0) * (m0))./ (k_0.* r* q_0)  - besselj(m0 + 1,   k_0.* r* q_0));
          
          Ephi_q_out_plus0  = @(r)  A0.* ((-p_0.* (m0./r)).* Ez_q0(r) + 1i * dHz_q0(r));
          Hphi_q_out_plus0  = @(r)  A0.* (-1i * dEz_q0(r) - (p_0.* (m0./r)).* Hz_q0(r));
          %     Hrho_q_out_plus0  = @(r)  A0.* (((m0./r)).* Ez_q0(r) - 1i * p_0 * dHz_q0(r));
          %     Erho_q_out_plus0  = @(r)  A0.* (-1i * p_0 * dEz_q0(r) - (m0./r).* Hz_q0(r));
          
          
          %     S_plus = integral(@(r) (Erho_q_out_plus0(r).* Hphi_q_out_minus(r) - Ephi_q_out_plus0(r).* Hrho_q_out_minus(r)-...
          %         Erho_q_out_minus(r).* Hphi_q_out_plus0(r) + Ephi_q_out_minus(r).* Hrho_q_out_plus0(r)).* r, a_0, 10 * a_0)
          S_plus = ( -Ez_q(a_0).* Hphi_q_out_plus0(a_0) + Ephi_q_out_minus(a_0).* Hz_q0(a_0) +...
              Ez_q0(a_0).* Hphi_q_out_minus(a_0) - Ephi_q_out_plus0(a_0).*  Hz_q(a_0)).* (a_0./(1i*k_0*(p + p_0)));
          %%%%%%%%%%% end test %%%%%%%%%%%%%%%%%%%%
          
          %     Ephi_q_plus0(a_0)
          % Ephi_q_out_plus0(a_0)
          % 
          %     Hphi_q_plus0(a_0)
          % Hphi_q_out_plus0(a_0)
          
          
          
          
          m = -m;
          p = -p;          
          
% %           m_old = m0;
          m = abs(m);
          %%% analytical formula
          [B_1_mn,B_2_mn,Cm2_mn, Dm2_mn] = coefficientsOfContinuousSpectrum(typeOfCylinder, k_0, k_0, a_0, EE1,0,0, MU1, EE, MU, -m, (-p), 0, 1);
          [B_1_pl,B_2_pl,Cm2_pl, Dm2_pl] = coefficientsOfContinuousSpectrum(typeOfCylinder, k_0, k_0, a_0, EE1,0,0, MU1, EE, MU, m, (p), 0, 1);
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
          
          N_n = (I_into - I_out);
          
          
          %     aplus2 = (a_0 * exp(-1i*(p+p_0)*k_0*z)./(1i*k_0*(p+p_0))).*...
          %         (S_plus - S_minus)./ N_n;
          
          aplus2 = (S_plus + S_minus)./ N_n;
          
      case 'Gyrotropic'
          [B_1, B_2, Cm2, Dm2]  = coefficientsOfContinuousSpectrum(typeOfCylinder, k_0, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU,  m,   p, q,  0, 1);
          
          
          %%%% inner integral
%           q1 = sqrt(EE1*MU1 -  p.^2);
%           q10 = sqrt(EE1*MU1 -  p_0.^2);
%           A_1 = 1./ (k_0 * (EE1*MU1 - EE*MU * p.^2));
          %     A_10 = 1./ (k_0 * (EE1*MU1 - EE*MU * p_0.^2)); 
          %     A = 1./ (k_0 * (1 - p.^2)); 
          
%           Ez_q   = @(r) B_1.* besselj(m, k_0.* r* q1);
%           Hz_q   = @(r) B_2.* besselj(m, k_0.* r* q1);
%           dEz_q  = @(r) B_1.* k_0.* q1 *((besselj(m, k_0.* r* q1) * m)./ (k_0.* r* q1)  - besselj(m + 1, k_0.* r* q1));
%           dHz_q  = @(r) B_2.* k_0.* q1 *((besselj(m, k_0.* r* q1) * m)./ (k_0.* r* q1)  - besselj(m + 1, k_0.* r* q1));
          
          
          [q1, q2, n1, n2, alp1, alp2, bet1, bet2] = term_of_gyrotropic_waveguide(EE1, GG1, HH1, p);
          
          JM1    = @(r) besselj(m+1, (q1.* r).* k_0);
          JM2    = @(r) besselj(m+1, (q2.* r).* k_0);
          Jm1    = @(r) besselj(m, (q1.* r).* k_0);
          Jm2    = @(r) besselj(m, (q2.* r).* k_0);
          Jm1_Q1 = @(r) Jm1(r)./ ((q1.* r).* k_0);
          Jm2_Q2 = @(r) Jm2(r)./ ((q2.* r).* k_0);
            
          Ephi_q_minus  = @(r)  B_1 * 1i * (JM1(r) + (alp1.* m).* Jm1_Q1(r))+...
                                B_2 * 1i * (JM2(r) + (alp2.* m).* Jm2_Q2(r));
          Hrho_q_minus  = @(r) - 1i * B_1 * (p*JM1(r) - n1.* (bet1 * m).* Jm1_Q1(r))+...
                               - 1i * B_2 * (p*JM2(r) - n2.* (bet2 * m).* Jm2_Q2(r));
          Hphi_q_minus  = @(r) - B_1 * n1.* (JM1(r) - (bet1 * m).* Jm1_Q1(r))+...
                               - B_2 * n2.* (JM2(r) - (bet2 * m).* Jm2_Q2(r));
          Erho_q_minus  = @(r) - B_1 * ((n1*p + GG1)/EE1 * JM1(r) + (alp1.* m).* Jm1_Q1(r))+...
                               - B_2 * ((n2*p + GG1)/EE1 * JM2(r) + (alp2.* m).* Jm2_Q2(r));
          
          %%%%% Bessel Beam
          A = 1./ (k_0 * (1 - p.^2));
          A0 = 1./ (k_0 * (1 - p_0.^2));
%           
          Ez_q0   = @(r) AE_0.* besselj(m0,  k_0.* r* q_0);
          Hz_q0   = @(r) AH_0.* besselj(m0,  k_0.* r* q_0);
          dEz_q0  = @(r) AE_0.* k_0.* (q_0) *((besselj(m0,   k_0.* r* q_0) * (m0))./ (k_0.* r* q_0)  - besselj(m0 + 1,   k_0.* r* q_0));
          dHz_q0  = @(r) AH_0.* k_0.* (q_0) *((besselj(m0,   k_0.* r* q_0) * (m0))./ (k_0.* r* q_0)  - besselj(m0 + 1,   k_0.* r* q_0));
          
          Ephi_q_plus0  = @(r)  A0.* ((-p_0.* (m0./r)).* Ez_q0(r) + 1i * dHz_q0(r));
          Hphi_q_plus0  = @(r)  A0.* (-1i * dEz_q0(r) - (p_0.* (m0./r)).* Hz_q0(r));
          Hrho_q_plus0  = @(r)  A0.* (((m0./r)).* Ez_q0(r) - 1i * p_0 * dHz_q0(r));
          Erho_q_plus0  = @(r)  A0.* (-1i * p_0 * dEz_q0(r) - (m0./r).* Hz_q0(r));
          
          
          S_minus = integral(@(r) (Erho_q_plus0(r).* Hphi_q_minus(r) - Ephi_q_plus0(r).* Hrho_q_minus(r)-...
              Erho_q_minus(r).* Hphi_q_plus0(r) + Ephi_q_minus(r).* Hrho_q_plus0(r)).* r, 0, a_0);
 
          
          %%%% outer integral   
          Ez_q   = @(r) Cm2.* q.* besselh(m,2, k_0.* r* q);
          Hz_q   = @(r) Dm2.* q.* besselh(m,2, k_0.* r* q);
          dEz_q  = @(r) Cm2.* q.* k_0.* q.* ((besselh(m,2, k_0.* r* q) * m)./ (k_0.* r* q)  - besselh(m + 1,2, k_0.* r* q));
          dHz_q  = @(r) Dm2.* q.* k_0.* q.* ((besselh(m,2, k_0.* r* q) * m)./ (k_0.* r* q)  - besselh(m + 1,2, k_0.* r* q));
          
          Ephi_q_out_minus  = @(r)  A.* ((-p.* (m./r)).* Ez_q(r) + 1i * dHz_q(r));
          %     Hrho_q_out_minus  = @(r)  A.* (((m./r)).* Ez_q(r) - 1i * p * dHz_q(r));
          Hphi_q_out_minus  = @(r)  A.* (-1i * dEz_q(r) - (p.* (m./r)).* Hz_q(r));
          %     Erho_q_out_minus  = @(r)  A.* (-1i * p * dEz_q(r) - (m./r).* Hz_q(r));
          
          
          %%%%% the Bessel Beam
          Ez_q0   = @(r) AE_0.* besselj(m0,  k_0.* r* q_0);
          Hz_q0   = @(r) AH_0.* besselj(m0,  k_0.* r* q_0);
          dEz_q0  = @(r) AE_0.* k_0.* (q_0) *((besselj(m0,   k_0.* r* q_0) * (m0))./ (k_0.* r* q_0)  - besselj(m0 + 1,   k_0.* r* q_0));
          dHz_q0  = @(r) AH_0.* k_0.* (q_0) *((besselj(m0,   k_0.* r* q_0) * (m0))./ (k_0.* r* q_0)  - besselj(m0 + 1,   k_0.* r* q_0));
          
          Ephi_q_out_plus0  = @(r)  A0.* ((-p_0.* (m0./r)).* Ez_q0(r) + 1i * dHz_q0(r));
          Hphi_q_out_plus0  = @(r)  A0.* (-1i * dEz_q0(r) - (p_0.* (m0./r)).* Hz_q0(r));
          %     Hrho_q_out_plus0  = @(r)  A0.* (((m0./r)).* Ez_q0(r) - 1i * p_0 * dHz_q0(r));
          %     Erho_q_out_plus0  = @(r)  A0.* (-1i * p_0 * dEz_q0(r) - (m0./r).* Hz_q0(r));
          
          
          %     S_plus = integral(@(r) (Erho_q_out_plus0(r).* Hphi_q_out_minus(r) - Ephi_q_out_plus0(r).* Hrho_q_out_minus(r)-...
          %         Erho_q_out_minus(r).* Hphi_q_out_plus0(r) + Ephi_q_out_minus(r).* Hrho_q_out_plus0(r)).* r, a_0, 10 * a_0)
          S_plus = ( -Ez_q(a_0).* Hphi_q_out_plus0(a_0) + Ephi_q_out_minus(a_0).* Hz_q0(a_0) +...
              Ez_q0(a_0).* Hphi_q_out_minus(a_0) - Ephi_q_out_plus0(a_0).*  Hz_q(a_0)).* (a_0./(1i*k_0*(p + p_0)));
          %%%%%%%%%%% end test %%%%%%%%%%%%%%%%%%%%

         m = -m;
          p = -p;
          GG1 = -GG1;
%           JM0    = besselj(m+1, (q_0.* a_0).* k_0);                            
%           JM1    = besselj(m+1, (q1.* a_0).* k_0);
%           JM2    = besselj(m+1, (q2.* a_0).* k_0);
%           JMM0   = besselj(m+2, (q_0.* a_0).* k_0);
%           JMM1   = besselj(m+2, (q1.* a_0).* k_0);
%           JMM2   = besselj(m+2, (q2.* a_0).* k_0);
%           Jm0    = besselj(m, (q_0.* a_0).* k_0);
%           Jm1    = besselj(m, (q1.* a_0).* k_0);
%           Jm2    = besselj(m, (q2.* a_0).* k_0);
%           
%           M10 =(a_0./ ((k_0*q1).^2 - (k_0*q_0).^2)).*...
%               (k_0*q1.* JMM1.* JM0 - k_0*q_0.* JMM0.* JM1);
%           M20 =(a_0./ ((k_0*q2).^2 - (k_0*q_0).^2)).*...
%               (k_0*q2.* JMM2.* JM0 - k_0*q_0.* JMM0.* JM2);
%           
%           A_0 = 1./ (k_0 * (1 - p_0.^2));
%           A_1 = 1./ (k_0 * (GG1^2 - (p.^2 - EE1).^2));
%           
%           I_into_1 = (p.* (EE1 - p.^2) - p_0* (GG1^2 - EE1.* (EE1 - p.^2)))*...
%               (1i/HH1)* AE_0.*...
%               ((B_1.* n1.* q1).* (m.* Jm1.* Jm0 + k_0^2 * q1.* q_0.* M10)+...
%                (B_2.* n2.* q2).* (m.* Jm2.* Jm0 + k_0^2 * q2.* q_0.* M20));
%           
%           I_into_2 = ((EE1 - p.^2).* (-p - p_0))*...
%               AH_0.*...
%               ((-B_1.* q1).* (m.* Jm1.* Jm0 + k_0^2 * q1.* q_0.* M10)+...
%                (-B_2.* q2).* (m.* Jm2.* Jm0 + k_0^2 * q2.* q_0.* M20));
%            
%           I_into_3 = m * (p.* p_0 + 1).* p.* GG1.*...
%                   ((1i./HH1).* (B_1.* n1.* q1.* Jm1 + B_2.* n2.* q2.* Jm2)).*...
%                                (AE_0.* Jm0);
%            
%           I_into_4 = GG1 * (-p - p_0).*...
%               m.* (-B_1.* q1.* Jm1 - B_2.* q2.* Jm2).* (AH_0.* Jm0);
%           
%           I_into_5 = 1i * (-p.* p_0 - 1).* (EE1 - p.^2).*...
%               m.* (-B_1.* q1.* Jm1 - B_2.* q2.* Jm2).* (AE_0.* Jm0);
%           
%           I_into_6 = 1i * (-p.* p_0.* (EE1 - p.^2) + (GG1^2 - EE1.* (EE1 - p.^2))).*...
%               m.* ((1i./HH1).* (B_1.* n1.* q1.* Jm1 + B_2.* n2.* q2.* Jm2)).* (AH_0.* Jm0);
% 
%           I_into_7 = 1i * GG1 * (-p.*p_0 - 1).*...
%                AE_0.*...
%               ((-B_1.* q1).* (m.* Jm1.* Jm0 + k_0^2 * q1.* q_0.* M10)+...
%                (-B_2.* q2).* (m.* Jm2.* Jm0 + k_0^2 * q2.* q_0.* M20));
%            
%            I_into_8 = 1i * GG1 * p.* (-p_0 - p).*...
%               (1i/HH1)* AH_0.*...
%               ((B_1.* n1.* q1).* (m.* Jm1.* Jm0 + k_0^2 * q1.* q_0.* M10)+...
%                (B_2.* n2.* q2).* (m.* Jm2.* Jm0 + k_0^2 * q2.* q_0.* M20));
%           
%           S_minus = (-1)^m * (A_1.* A_0).* (I_into_1 + I_into_2 + I_into_3 +...
%                                             I_into_4 + I_into_5 + I_into_6 +...
%                                             I_into_7 + I_into_8)
          %%% analytical formula
%           [B_1_mn,B_2_mn,Cm2_mn, Dm2_mn] = coefficientsOfContinuousSpectrum(typeOfCylinder, k_0, k_0, a_0, EE1,-GG1, HH1, MU1, EE, MU, -m, (-p), 0, 1);
%           [B_1_pl,B_2_pl,Cm2_pl, Dm2_pl] = coefficientsOfContinuousSpectrum(typeOfCylinder, k_0, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, (p), 0, 1);
%           %%% norm
%           q1 = sqrt(EE1*MU1/(EE*MU) -  p.^2);
%           A_1 = 1./ (k_0 * (EE1*MU1 - EE*MU * p.^2));       
%           A = 1./ (k_0 * (1 - p.^2)); 
%           
%           G_Jm_a0 = - 1/2 * (besselj(0, k_0*q1*a_0)).^2;
%           G_Jm_0  = - 1/2 * (besselj(0, k_0*q1*0)).^2;    
%           G_Hm2_a0 = - 1/2 * (besselh(0,2, k_0*q*a_0)).^2;
%           for n=0:(m-1)
%               G_Jm_a0  =  G_Jm_a0  - (besselj(m-n, k_0*q1*a_0)).^2;
%               G_Hm2_a0 =  G_Hm2_a0 - (besselh(m-n, 2, k_0*q*a_0)).^2;
%           end
%           
%           G_Jm_1__a0 = - 1/2 * (besselj(0, k_0*q1*a_0)).^2;
%           G_Jm_1__0  = - 1/2 * (besselj(0, k_0*q1*0)).^2;
%           G_Hm2_1__a0 = - 1/2 * (besselh(0,2, k_0*q*a_0)).^2;    
%           for n=0:(m-2)
%               G_Jm_1__a0  =  G_Jm_1__a0  - (besselj(m-1-n, k_0*q1*a_0)).^2;
%               G_Hm2_1__a0 =  G_Hm2_1__a0 - (besselh(m-1-n, 2, k_0*q*a_0)).^2;
%           end
%           
%           F_Jm_a0 = ((k_0*q1*a_0).^2 / 2).* ((besselj(m+1, k_0*q1*a_0)).^2  - besselj(m,   k_0*q1*a_0).* besselj(m+2,  k_0*q1*a_0));
%           F_Hm2_a0 =((k_0*q *a_0).^2 / 2).* ((besselh(m+1,2, k_0*q*a_0)).^2 - besselh(m,2, k_0*q *a_0).* besselh(m+2,2, k_0*q *a_0));
%           
%           I_into = (-1)^m * (A_1.^2).* (-2 * p)* (EE1 * B_1_mn.* B_1_pl - B_2_mn.* B_2_pl).*...
%               (-m * (G_Jm_a0 - G_Jm_0) + m * (G_Jm_1__a0 - G_Jm_1__0) + F_Jm_a0);
%           I_into = I_into + (-1)^m * (A_1.^2) * m * 1i * (p.^2 + EE1*MU1) *...
%               ((B_1_mn.* B_2_pl + B_2_mn.* B_1_pl).* (besselj(m, k_0*q1*a_0)).^2);%%% добавление полного дифференциала
%           
%           
%           I_out  =  (-1)^m * (A.^2).* (-2 * p)*(Cm2_mn.* Cm2_pl - Dm2_mn.* Dm2_pl).*...
%               (-m * G_Hm2_a0 + m * G_Hm2_1__a0 + F_Hm2_a0);
%           I_out  = I_out + (-1)^m * (A.^2) * m * 1i * (p.^2 + 1) *...
%               ((Cm2_mn.* Dm2_pl + Cm2_pl.* Dm2_mn).* (besselh(m,2, k_0*q*a_0)).^2);%%% добавление полного дифференциала
%           
%           N_n = (I_into - q.^2.* I_out);
%           
%           
%           
                
          
          
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
   
   
   
   
   
   
   
   
%    [q1, q2, n1, n2, alp1, alp2, bet1, bet2] = term_of_gyrotropic_waveguide(EE1, GG1, HH1, p);
%    [B_1,B_2,Cm2, Dm2] = coefficientsOfContinuousSpectrum('Gyrotropic', k_0, k_0, a_0, EE1, GG1,HH1, MU1, EE, MU,  m, (p), 0, 1);
%             Q1 = @(r)  (q1.* r).* k_0;
%             Q2 = @(r)  (q2.* r).* k_0;  
%             
%             JM1    =@(r)  besselj(m+1, Q1(r));
%             JM2    =@(r)  besselj(m+1, Q2(r));
%             Jm1    =@(r)  besselj(m, Q1(r));
%             Jm2    =@(r)  besselj(m, Q2(r));
%             Jm1_Q1 =@(r)  Jm1(r)./ Q1(r);
%             Jm2_Q2 =@(r)  Jm2(r)./ Q2(r);
%             
%             Ez_q1   =@(r)  (((1i./HH1).* n1).* q1).* Jm1(r);
%             Ez_q2   =@(r)  (((1i./HH1).* n2).* q2).* Jm2(r);
%             Hz_q1   =@(r)  - q1.* Jm1(r);
%             Hz_q2   =@(r)  - q2.* Jm2(r);            
%             dEz_q1   =@(r)  (((1i./HH1).* n1).* q1).* (m.* Jm1_Q1(r) - JM1(r)).* q1.* k_0;
%             dEz_q2   =@(r)  (((1i./HH1).* n2).* q2).* (m.* Jm2_Q2(r) - JM2(r)).* q2.* k_0;
%             dHz_q1   =@(r)  - q1.* (m.* Jm1_Q1(r) - JM1(r)).* q1.* k_0;
%             dHz_q2   =@(r)  - q2.* (m.* Jm2_Q2(r) - JM2(r)).* q2.* k_0;
%             A_1 = 1./ (k_0 * (GG1^2 - (p.^2 - EE1).^2));
%    
%             Ephi = @(r) A_1*(p.*(EE1-p.^2).*m./r.* (B_1.* Ez_q1(r)+ B_2.* Ez_q2(r))+...
%                                       p.*GG1.* (B_1.* dEz_q1(r)+ B_2.* dEz_q2(r))+...
%                                    -1i*GG1.*m./r.* (B_1.* Hz_q1(r)+ B_2.* Hz_q2(r))+...
%                                  -1i.*(EE1-p.^2).* (B_1.* dHz_q1(r)+ B_2.* dHz_q2(r)));
%                              
%             Erho = @(r) A_1*(1i*p.*GG1.*m./r.* (B_1.* Ez_q1(r)+  B_2.* Ez_q2(r))+...
%                             1i*p.*(EE1-p.^2).* (B_1.* dEz_q1(r)+ B_2.* dEz_q2(r))+...
%                                    (EE1-p.^2).*m./r.* (B_1.* Hz_q1(r)+ B_2.* Hz_q2(r))+...
%                                                 GG1.* (B_1.* dHz_q1(r)+ B_2.* dHz_q2(r)));
%             
%             Hrho = @(r) A_1*((GG1^2 - EE1*(EE1 - p.^2)).*m./r.* (B_1.* Ez_q1(r)+ B_2.* Ez_q2(r))+...
%                                       -p^2.*GG1.* (B_1.* dEz_q1(r)+ B_2.* dEz_q2(r))+...
%                                     1i*p*GG1.*m./r.* (B_1.* Hz_q1(r)+ B_2.* Hz_q2(r))+...
%                                   1i.*p*(EE1-p.^2).* (B_1.* dHz_q1(r)+ B_2.* dHz_q2(r)));
%                              
%             Hphi = @(r) A_1*(1i*p^2.*GG1.*m./r.* (B_1.* Ez_q1(r)+  B_2.* Ez_q2(r))+...
%                             -1i*(GG1^2 - EE1*(EE1 - p.^2)).* (B_1.* dEz_q1(r)+ B_2.* dEz_q2(r))+...
%                                   p*(EE1-p.^2).*m./r.* (B_1.* Hz_q1(r)+ B_2.* Hz_q2(r))+...
%                                               p*GG1.* (B_1.* dHz_q1(r)+ B_2.* dHz_q2(r)));
%                                             
% I_into = 2*integral(@(r) (-Erho(r).* Hphi(r) - Ephi(r).* Hrho(r)).* r, 0, a_0);
   
   
   
 m_Old = m;
 m = abs(m);

%%% norm
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
    
    N_n = (I_into - I_out);
          
          
          %     aplus2 = (a_0 * exp(-1i*(p+p_0)*k_0*z)./(1i*k_0*(p+p_0))).*...
          %         (S_plus - S_minus)./ N_n;
          
          aplus2 = (S_plus + S_minus)./ N_n;
          
  end





