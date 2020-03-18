clear all
tic
systemParameters
q_0 = -1.3836 * 1i;

% p_0 = 1/sqrt(2);
% q_0 = sqrt(1-p_0^2);
% q_0 = q_0.* (2*(imag(q_0) <= 0)-1);

%%%%%% paramenter of cylinder for perfect conductivity core
% a_0 = 2.4048./ (k_0 * q_0);

m=1
p_0 = sqrt(1-q_0^2);
p_0 = real(p_0) - 1i * abs(imag(p_0));

p = [ 
        2.2196;
        ]; %%%k_0 = pi / a_0/ 2;
p = [  2.66233]; %%% GG1 = 1; k_0 = pi / a_0/ 2;

p = [ 3.176710633456818];

p = [2.84700547735529];

p_n__of_descreteMode_of_gyrotropicCyl;
p = p_n(1);


% p_n = 3.1622566;  %%%k_0 = pi / a_0/ 10;
q = sqrt(1-p.^2);
q = real(q) - 1i * abs(imag(q));

% GG1 = 100
% HH1 = -10;
% GG1=-GG1
% p = -p;
% m = -m;

% p=-p;


[B_1_mn,B_2_mn,Cm2_mn, Dm2_mn] = coefficientsOfContinuousSpectrum('Gyrotropic', k_0, k_0, a_0, EE1,-GG1,HH1, MU1, EE, MU, -m, (-p), q, 0, 1);
[B_1_pl,B_2_pl,Cm2_pl, Dm2_pl] = coefficientsOfContinuousSpectrum('Gyrotropic', k_0, k_0, a_0, EE1, GG1,HH1, MU1, EE, MU,  m, ( p), q, 0, 1);


%%%% inner integral
          [q1, q2, n1, n2, alp1, alp2, bet1, bet2] = term_of_gyrotropic_waveguide(EE1, -GG1, HH1, -p);

          JM1    = @(r) besselj(m+1, (q1.* r).* k_0);
          JM2    = @(r) besselj(m+1, (q2.* r).* k_0);
          Jm1    = @(r) besselj(m, (q1.* r).* k_0);
          Jm2    = @(r) besselj(m, (q2.* r).* k_0);
          Jm1_Q1 = @(r) Jm1(r)./ ((q1.* r).* k_0);
          Jm2_Q2 = @(r) Jm2(r)./ ((q2.* r).* k_0);
          
          JM1_mn    = @(r) besselj(-m+1, (q1.* r).* k_0);
          JM2_mn    = @(r) besselj(-m+1, (q2.* r).* k_0);
          Jm1_mn    = @(r) besselj(-m, (q1.* r).* k_0);
          Jm2_mn    = @(r) besselj(-m, (q2.* r).* k_0);
          Jm1_Q1_mn = @(r) Jm1_mn(r)./ ((q1.* r).* k_0);
          Jm2_Q2_mn = @(r) Jm2_mn(r)./ ((q2.* r).* k_0);

          
          Ez_q_minus = @(r) B_1_mn.*  (((1i./HH1).* n1).* q1).* Jm1_mn(r) +...
                            B_2_mn.*  (((1i./HH1).* n2).* q2).* Jm2_mn(r);
          Hz_q_minus  = @(r) - B_1_mn.*  q1.* Jm1_mn(r) +...
                             - B_2_mn.*  q2.* Jm2_mn(r);
          Ephi_q_minus  = @(r)  B_1_mn * 1i * (JM1_mn(r) + (alp1.* (-m)).* Jm1_Q1_mn(r))+...
                                B_2_mn * 1i * (JM2_mn(r) + (alp2.* (-m)).* Jm2_Q2_mn(r));
          Hrho_q_minus  = @(r) - 1i * B_1_mn * (-p.* JM1_mn(r) - n1.* (bet1 * (-m)).* Jm1_Q1_mn(r))+...
                               - 1i * B_2_mn * (-p.* JM2_mn(r) - n2.* (bet2 * (-m)).* Jm2_Q2_mn(r));
          Hphi_q_minus  = @(r) - B_1_mn * n1.* (JM1_mn(r) - (bet1 * (-m)).* Jm1_Q1_mn(r))+...
                               - B_2_mn * n2.* (JM2_mn(r) - (bet2 * (-m)).* Jm2_Q2_mn(r));
          Erho_q_minus  = @(r) - B_1_mn * ((-n1*p - GG1)/EE1 * JM1_mn(r) + (alp1.* (-m)).* Jm1_Q1_mn(r))+...
                               - B_2_mn * ((-n2*p - GG1)/EE1 * JM2_mn(r) + (alp2.* (-m)).* Jm2_Q2_mn(r));

          [q1, q2, n1_pl, n2_pl, alp1_pl, alp2_pl, bet1_pl, bet2_pl] = term_of_gyrotropic_waveguide(EE1, GG1, HH1, p);                 

          Ez_q_plus  = @(r) B_1_pl.*  (((1i./HH1).* n1_pl).* q1).* Jm1(r) +...
                            B_2_pl.*  (((1i./HH1).* n2_pl).* q2).* Jm2(r);
          Hz_q_plus  = @(r) - B_1_pl.*  q1.* Jm1(r) +...
                            - B_2_pl.*  q2.* Jm2(r);
          Ephi_q_plus  = @(r)  B_1_pl * 1i * (JM1(r) + (alp1_pl.* m).* Jm1_Q1(r))+...
                               B_2_pl * 1i * (JM2(r) + (alp2_pl.* m).* Jm2_Q2(r));
          Hrho_q_plus  = @(r) - 1i * B_1_pl * (p*JM1(r) - n1_pl.* (bet1_pl * m).* Jm1_Q1(r))+...
                              - 1i * B_2_pl * (p*JM2(r) - n2_pl.* (bet2_pl * m).* Jm2_Q2(r));
          Hphi_q_plus  = @(r) - B_1_pl * n1_pl.* (JM1(r) - (bet1_pl * m).* Jm1_Q1(r))+...
                              - B_2_pl * n2_pl.* (JM2(r) - (bet2_pl * m).* Jm2_Q2(r));
          Erho_q_plus  = @(r) - B_1_pl * ((n1_pl*(p) + GG1)/EE1 * JM1(r) + (alp1_pl.* m).* Jm1_Q1(r))+...
                              - B_2_pl * ((n2_pl*(p) + GG1)/EE1 * JM2(r) + (alp2_pl.* m).* Jm2_Q2(r));

    I_into_num = integral(@(r) (Erho_q_plus(r).* Hphi_q_minus(r) - Ephi_q_plus(r).* Hrho_q_minus(r)-...
                               Erho_q_minus(r).* Hphi_q_plus(r)  + Ephi_q_minus(r).* Hrho_q_plus(r)).* r, 0, a_0)
                           
                           
                           
                           
                           
          Ephi_q_minus  = @(r)  B_1_mn * 1i * (JM1_mn(r) + (alp1.* (-m)).* Jm1_Q1_mn(r))+...
                                B_2_mn * 1i * (JM2_mn(r) + (alp2.* (-m)).* Jm2_Q2_mn(r));
          Hrho_q_minus  = @(r) - 1i * B_1_mn * (-p.* JM1_mn(r) - n1.* (bet1 * (-m)).* Jm1_Q1_mn(r))+...
                               - 1i * B_2_mn * (-p.* JM2_mn(r) - n2.* (bet2 * (-m)).* Jm2_Q2_mn(r));
          Hphi_q_minus  = @(r) - B_1_mn * n1.* (JM1_mn(r) - (bet1 * (-m)).* Jm1_Q1_mn(r))+...
                               - B_2_mn * n2.* (JM2_mn(r) - (bet2 * (-m)).* Jm2_Q2_mn(r));
          Erho_q_minus  = @(r) - B_1_mn * ((-n1*p - GG1)/EE1 * JM1_mn(r) + (alp1.* (-m)).* Jm1_Q1_mn(r))+...
                               - B_2_mn * ((-n2*p - GG1)/EE1 * JM2_mn(r) + (alp2.* (-m)).* Jm2_Q2_mn(r));

          [q1, q2, n1_pl, n2_pl, alp1_pl, alp2_pl, bet1_pl, bet2_pl] = term_of_gyrotropic_waveguide(EE1, GG1, HH1, p);           

          Ephi_q_plus  = @(r)  B_1_pl * 1i * (JM1(r) + (alp1_pl.* m).* Jm1_Q1(r))+...
                               B_2_pl * 1i * (JM2(r) + (alp2_pl.* m).* Jm2_Q2(r));
          Hrho_q_plus  = @(r) - 1i * B_1_pl * (p*JM1(r) - n1_pl.* (bet1_pl * m).* Jm1_Q1(r))+...
                              - 1i * B_2_pl * (p*JM2(r) - n2_pl.* (bet2_pl * m).* Jm2_Q2(r));
          Hphi_q_plus  = @(r) - B_1_pl * n1_pl.* (JM1(r) - (bet1_pl * m).* Jm1_Q1(r))+...
                              - B_2_pl * n2_pl.* (JM2(r) - (bet2_pl * m).* Jm2_Q2(r));
          Erho_q_plus  = @(r) - B_1_pl * ((n1_pl*(p) + GG1)/EE1 * JM1(r) + (alp1_pl.* m).* Jm1_Q1(r))+...
                              - B_2_pl * ((n2_pl*(p) + GG1)/EE1 * JM2(r) + (alp2_pl.* m).* Jm2_Q2(r));

    I_into_num2 = integral(@(r) (Erho_q_plus(r).* Hphi_q_minus(r) - Ephi_q_plus(r).* Hrho_q_minus(r)-...
                               Erho_q_minus(r).* Hphi_q_plus(r)  + Ephi_q_minus(r).* Hrho_q_plus(r)).* r, 0, a_0)
    I_into_num22 = 2*integral(@(r) (-Erho_q_plus(r).* Hphi_q_plus(r) - Ephi_q_plus(r).* Hrho_q_plus(r)).* r, 0, a_0)                           
                           
    

  %%%% from Gauss theorem
p_0 = 2.219599;
p_0 = [  2.6623299]; %%% GG1 = -1; k_0 = pi / a_0/ 2;
p_0 = [2.84700548];

p_0 = p + 1e-8;
q_0 = sqrt(1-p_0.^2);
q_0 = real(q_0) - 1i * abs(imag(q_0));
    
    [q1, q2, n1_pl, n2_pl, alp1_pl, alp2_pl, bet1_pl, bet2_pl] = term_of_gyrotropic_waveguide(EE1, GG1, HH1, p_0);                 
     JM1    = @(r) besselj(m+1, (q1.* r).* k_0);
     JM2    = @(r) besselj(m+1, (q2.* r).* k_0);
     Jm1    = @(r) besselj(m, (q1.* r).* k_0);
     Jm2    = @(r) besselj(m, (q2.* r).* k_0);
     Jm1_Q1 = @(r) Jm1(r)./ ((q1.* r).* k_0);
     Jm2_Q2 = @(r) Jm2(r)./ ((q2.* r).* k_0);
  
     [B_1_pl,B_2_pl,Cm2_pl, Dm2_pl] = coefficientsOfContinuousSpectrum('Gyrotropic', k_0, k_0, a_0, EE1, GG1,HH1, MU1, EE, MU, m, (p_0), q_0, 0, 1);
          Ez_q_0  = @(r) B_1_pl.*  (((1i./HH1).* n1_pl).* q1).* Jm1(r) +...
                         B_2_pl.*  (((1i./HH1).* n2_pl).* q2).* Jm2(r);
          Hz_q_0  = @(r) - B_1_pl.*  q1.* Jm1(r) +...
                         - B_2_pl.*  q2.* Jm2(r);
    Ephi_q_plus0  = @(r)  B_1_pl * 1i * (JM1(r) + (alp1_pl.* m).* Jm1_Q1(r))+...
                          B_2_pl * 1i * (JM2(r) + (alp2_pl.* m).* Jm2_Q2(r));
    Hrho_q_plus0  = @(r) - 1i * B_1_pl * (p_0*JM1(r) - n1_pl.* (bet1_pl * m).* Jm1_Q1(r))+...
                         - 1i * B_2_pl * (p_0*JM2(r) - n2_pl.* (bet2_pl * m).* Jm2_Q2(r));
    Hphi_q_plus0  = @(r) - B_1_pl * n1_pl.* (JM1(r) - (bet1_pl * m).* Jm1_Q1(r))+...
                         - B_2_pl * n2_pl.* (JM2(r) - (bet2_pl * m).* Jm2_Q2(r));
    Erho_q_plus0  = @(r) - B_1_pl * ((n1_pl*(p_0) + GG1)/EE1 * JM1(r) + (alp1_pl.* m).* Jm1_Q1(r))+...
                         - B_2_pl * ((n2_pl*(p_0) + GG1)/EE1 * JM2(r) + (alp2_pl.* m).* Jm2_Q2(r));

    
    I_into_jump = (-Ez_q_minus(a_0).*  Hphi_q_plus0(a_0) +...
                    Ephi_q_minus(a_0).* Hz_q_0(a_0) +...
                    Ez_q_0(a_0).* Hphi_q_minus(a_0) -...
                    Ephi_q_plus0(a_0).* Hz_q_minus(a_0));
            
    I_into_jump = I_into_jump.* (a_0./(1i*k_0*(-p + p_0)))  
    

    
%%% analytical formula
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
   M12 =(a_0./ ((k_0*q1)^2 - (k_0*q2)^2)).*...
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
                                     I_into_4 + I_into_5 + I_into_6)

   
% Ez11 = (-1)^m * (A_1.^2).* 2.* (I_into_1 + I_into_4)
% Hz11 = (-1)^m * (A_1.^2).* 2.* (I_into_2 + I_into_5)
% EHz11 = (-1)^m * (A_1.^2).* 2.* (I_into_3)

   
   [q1, q2, n1, n2, alp1, alp2, bet1, bet2] = term_of_gyrotropic_waveguide(EE1, GG1, HH1, p);
   [B_1,B_2,Cm2, Dm2] = coefficientsOfContinuousSpectrum('Gyrotropic', k_0, k_0, a_0, EE1, GG1,HH1, MU1, EE, MU,  m, (p), q, 0, 1);
            Q1 = @(r)  (q1.* r).* k_0;
            Q2 = @(r)  (q2.* r).* k_0;  
            
            JM1    =@(r)  besselj(m+1, Q1(r));
            JM2    =@(r)  besselj(m+1, Q2(r));
            Jm1    =@(r)  besselj(m, Q1(r));
            Jm2    =@(r)  besselj(m, Q2(r));
            Jm1_Q1 =@(r)  Jm1(r)./ Q1(r);
            Jm2_Q2 =@(r)  Jm2(r)./ Q2(r);
            
            Ez_q1   =@(r)  (((1i./HH1).* n1).* q1).* Jm1(r);
            Ez_q2   =@(r)  (((1i./HH1).* n2).* q2).* Jm2(r);
            Hz_q1   =@(r)  - q1.* Jm1(r);
            Hz_q2   =@(r)  - q2.* Jm2(r);            
            dEz_q1   =@(r)  (((1i./HH1).* n1).* q1).* (m.* Jm1_Q1(r) - JM1(r)).* q1.* k_0;
            dEz_q2   =@(r)  (((1i./HH1).* n2).* q2).* (m.* Jm2_Q2(r) - JM2(r)).* q2.* k_0;
            dHz_q1   =@(r)  - q1.* (m.* Jm1_Q1(r) - JM1(r)).* q1.* k_0;
            dHz_q2   =@(r)  - q2.* (m.* Jm2_Q2(r) - JM2(r)).* q2.* k_0;
            A_1 = 1./ (k_0 * (GG1^2 - (p.^2 - EE1).^2));
   
            Ephi = @(r) A_1*(p.*(EE1-p.^2).*m./r.* (B_1.* Ez_q1(r)+ B_2.* Ez_q2(r))+...
                                      p.*GG1.* (B_1.* dEz_q1(r)+ B_2.* dEz_q2(r))+...
                                   -1i*GG1.*m./r.* (B_1.* Hz_q1(r)+ B_2.* Hz_q2(r))+...
                                 -1i.*(EE1-p.^2).* (B_1.* dHz_q1(r)+ B_2.* dHz_q2(r)));
                             
            Erho = @(r) A_1*(1i*p.*GG1.*m./r.* (B_1.* Ez_q1(r)+  B_2.* Ez_q2(r))+...
                            1i*p.*(EE1-p.^2).* (B_1.* dEz_q1(r)+ B_2.* dEz_q2(r))+...
                                   (EE1-p.^2).*m./r.* (B_1.* Hz_q1(r)+ B_2.* Hz_q2(r))+...
                                                GG1.* (B_1.* dHz_q1(r)+ B_2.* dHz_q2(r)));
            
            Hrho = @(r) A_1*((GG1^2 - EE1*(EE1 - p.^2)).*m./r.* (B_1.* Ez_q1(r)+ B_2.* Ez_q2(r))+...
                                      -p^2.*GG1.* (B_1.* dEz_q1(r)+ B_2.* dEz_q2(r))+...
                                    1i*p*GG1.*m./r.* (B_1.* Hz_q1(r)+ B_2.* Hz_q2(r))+...
                                  1i.*p*(EE1-p.^2).* (B_1.* dHz_q1(r)+ B_2.* dHz_q2(r)));
                             
            Hphi = @(r) A_1*(1i*p^2.*GG1.*m./r.* (B_1.* Ez_q1(r)+  B_2.* Ez_q2(r))+...
                            -1i*(GG1^2 - EE1*(EE1 - p.^2)).* (B_1.* dEz_q1(r)+ B_2.* dEz_q2(r))+...
                                  p*(EE1-p.^2).*m./r.* (B_1.* Hz_q1(r)+ B_2.* Hz_q2(r))+...
                                              p*GG1.* (B_1.* dHz_q1(r)+ B_2.* dHz_q2(r)));
                                            
I_into_num23 = 2*integral(@(r) (-Erho(r).* Hphi(r) - Ephi(r).* Hrho(r)).* r, 0, a_0)

% I_into_num23-Ez11-Hz11
% 
% (-1)^m * (A_1.^2).* 2.* (1i* (-(EE1 - p.^2).* (GG1^2 - EE1.* (EE1 - p.^2)) +...
%                      p.^2 * (EE1 - p.^2).^2 +...
%                      2 * p.^2* GG1^2).*...
%        m.* (B_1.* Ez_q1(a_0) +  B_2.* Ez_q2(a_0)).*...
%                       (B_1.* Hz_q1(a_0) +  B_2.* Hz_q2(a_0)))
% 
%
%    I_into_3 = (-1i * (EE1 - p.^2).* (GG1^2 - EE1.* (EE1 - p.^2)) +...
%                 1i * p.^2 * (EE1 - p.^2).^2 + 2*1i* p.^2* GG1^2)*...
%        m * (B_1.* Ez_q1(a_0) + B_2.* Ez_q2(a_0)).*...
%            (B_1.* Hz_q1(a_0) + B_2.* Hz_q2(a_0))
%                   
%    I_into_4 = (p.*GG1.*(GG1^2 - EE1.* (EE1 - p.^2)) - p.^3 * GG1 * (EE1 - p.^2))*...
%        m * ((B_1.* Ez_q1(a_0) + B_2.* Ez_q2(a_0))).^2
%                      
%    I_into_5 = 2*(p.*GG1.* (EE1 - p.^2))*...
%        m * (B_1.* Hz_q1(a_0) + B_2.* Hz_q2(a_0)).^2
   
%    M12
%    F2
%    (k_0*q2).^2 * integral(@(r) besselj(m+1,k_0*q2*r).* besselj(m+1,k_0*q2*r).* r,0,a_0)
%        ((B_1.* q1).^2 * (m.* Jm1.^2 + F1)+...
%         (B_2.* q2).^2 * (m.* Jm2.^2 + F2)+...
%        2*(B_1.* q1).* (B_2.* q2).* (m.* Jm1.* Jm2 + k_0^2 * q1.* q2.* M12))
%    
%  integral(@(r) ((m./r.* ((B_1.* q1).* besselj(m,k_0*q1*r) + (B_2.* q2).* besselj(m,k_0*q2*r))).^2 + ...
%      (m./r.* ((B_1.* q1).* besselj(m,k_0*q1*r) + (B_2.* q2).* besselj(m,k_0*q2*r)) -...
%      ((B_1.* q1).* k_0*q1* besselj(m+1,k_0*q1*r) + (B_2.* q2).*k_0*q2* besselj(m+1,k_0*q2*r))).^2).* r,0,a_0)  
   



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
    
    I_out  =  (-1)^m * (A.^2).* (-2 * p)*(Cm2_mn.* Cm2_pl - Dm2_mn.* Dm2_pl).*...
        (-m * G_Hm2_a0 + m * G_Hm2_1__a0 + F_Hm2_a0);
    I_out  = I_out + (-1)^m * (A.^2) * m * 1i * (p.^2 + 1) *...
        ((Cm2_mn.* Dm2_pl + Cm2_pl.* Dm2_mn).* (besselh(m,2, k_0*q*a_0)).^2);
    
    I_out = q.^2 * I_out;
    
    N_n = (I_into - I_out);
    
    
    
    
    %%%% norm is calculated with using the new simplify formulas
    I_out_simplify = (-1)^(m+1) * (p * a_0^2/2 * (Dm2_mn^2 - Cm2_mn^2).*...
    ((besselh(m+1,2, k_0*q*a_0)).^2 - besselh(m,2, k_0*q *a_0).* besselh(m+2,2, k_0*q *a_0) +...
    2*m./(k_0*q*a_0)^2 * besselh(m,2, k_0*q *a_0)^2) +...
    1i* m * a_0^2./ (k_0*q*a_0)^2 * (p.^2 + 1) * Cm2_mn * Dm2_mn * besselh(m,2, k_0*q *a_0)^2)      *2;

  

    JM1    = besselj(m+1, (q1.* a_0).* k_0);
    JM2    = besselj(m+1, (q2.* a_0).* k_0);
    JMM1    = besselj(m+2, (q1.* a_0).* k_0);
    JMM2    = besselj(m+2, (q2.* a_0).* k_0);
    Jm1    = besselj(m, (q1.* a_0).* k_0);
    Jm2    = besselj(m, (q2.* a_0).* k_0);
    Q1 = q1.* a_0.* k_0;
    Q2 = q2.* a_0.* k_0;
       
    I1_1 = (a_0.^2 / 2).* (JM1.^2  - Jm1.*  JMM1);
    I1_2 = (a_0.^2 / 2).* (JM2.^2  - Jm2.*  JMM2);
    I_4 =  (a_0.^2./ (Q1^2 - Q2^2)).*...
                 (Q1.* JMM1.* JM2 - Q2.* JMM2.* JM1);
             
    I_in_simplify = (-1)^(m).*... 
                    (B_1^2.* ( I1_1.* (p + n1.* (n1.* p + GG1)./ EE1) - m./ (k_0* q1)^2 * alp1.* bet1.* n1.* Jm1.^2) +...
                     B_2^2.* ( I1_2.* (p + n2.* (n2.* p + GG1)./ EE1) - m./ (k_0* q2)^2 * alp2.* bet2.* n2.* Jm2.^2) +...
                     B_1.* B_2.* (I_4.* (2 * p * n1 * n2/ EE1 + 2 * p + GG1/EE1 * (n2 + n1)) +...
                     - m.* a_0^2/ (Q1.* Q2).* (alp1.* bet2.* n2 + alp2.* bet1.* n1).* Jm1.* Jm2)) * 2
    
    
    
%  p = [ 
%         2.2196;
%         ]   ;
%     q = sqrt(1-p.^2);
% q = real(q) - 1i * abs(imag(q));
% 
% p=-p;
% 
% [B_1_mn,B_2_mn,Cm2_mn, Dm2_mn] = coefficientsOfContinuousSpectrum('Isotropic',k_0, k_0, a_0, EE1,0,0, MU1, EE, MU, -m, (-p), 0, 1);
% [B_1_pl,B_2_pl,Cm2_pl, Dm2_pl] = coefficientsOfContinuousSpectrum('Isotropic',k_0, k_0, a_0, EE1,0,0, MU1, EE, MU, m, (p), 0, 1);
% 
% %%%% inner integral
%     q1 = sqrt(EE1*MU1/(EE*MU) -  p.^2);
%     A_1 = 1./ (k_0 * (EE1*MU1 - EE*MU * p.^2));       
%     A = 1./ (k_0 * (1 - p.^2)); 
%     
%     Ez_q_mn   = @(r) B_1_mn.* besselj(-m, k_0.* r* q1);
%     Hz_q_mn   = @(r) B_2_mn.* besselj(-m, k_0.* r* q1);
%     dEz_q_mn  = @(r) B_1_mn.* k_0.* q1 *((besselj(-m, k_0.* r* q1) * (-m))./ (k_0.* r* q1)  - besselj(-m + 1, k_0.* r* q1));
%     dHz_q_mn  = @(r) B_2_mn.* k_0.* q1 *((besselj(-m, k_0.* r* q1) * (-m))./ (k_0.* r* q1)  - besselj(-m + 1, k_0.* r* q1));
%     
%     Ephi_q_minus  = @(r)  A_1.* (( p.* (-m./r)).* Ez_q_mn(r) + 1i * dHz_q_mn(r));
%     Hrho_q_minus  = @(r)  A_1.* ((EE1 * (-m./r)).* Ez_q_mn(r) + 1i * p * dHz_q_mn(r));
%     Hphi_q_minus  = @(r)  A_1.* (-1i * EE1 * dEz_q_mn(r) + (p.* (-m./r)).* Hz_q_mn(r));
%     Erho_q_minus  = @(r)  A_1.* ( 1i * p * dEz_q_mn(r) - (-m./r).* Hz_q_mn(r));
%     
% 
%     Ez_q   = @(r) B_1_pl.* besselj(m, k_0.* r* q1);
%     Hz_q   = @(r) B_2_pl.* besselj(m, k_0.* r* q1);
%     dEz_q  = @(r) B_1_pl.* k_0.* q1 *((besselj(m, k_0.* r* q1) * m)./ (k_0.* r* q1)  - besselj(m + 1, k_0.* r* q1));
%     dHz_q  = @(r) B_2_pl.* k_0.* q1 *((besselj(m, k_0.* r* q1) * m)./ (k_0.* r* q1)  - besselj(m + 1, k_0.* r* q1));
%     
%     Ephi_q_plus  = @(r)  A_1.* ((-p.* (m./r)).* Ez_q(r) + 1i * dHz_q(r));
%     Hrho_q_plus  = @(r)  A_1.* ((EE1 * (m./r)).* Ez_q(r) - 1i * p * dHz_q(r));
%     Hphi_q_plus  = @(r)  A_1.* (-1i * EE1 * dEz_q(r) + (-p.* (m./r)).* Hz_q(r));
%     Erho_q_plus  = @(r)  A_1.* (-1i * p * dEz_q(r) - (m./r).* Hz_q(r));
% 
% 
%     I_into_num = integral(@(r) (Erho_q_plus(r).* Hphi_q_minus(r) - Ephi_q_plus(r).* Hrho_q_minus(r)-...
%                                Erho_q_minus(r).* Hphi_q_plus(r) + Ephi_q_minus(r).* Hrho_q_plus(r)).* r, 0, a_0)
%     
% %%%% from Gauss theorem
% p_0 = 2.21959;
% q_0 = sqrt(1-p_0^2);
% q_0 = q_0.* (2*(imag(q_0) <= 0)-1);
% 
% [B_1,B_2,Cm2, Dm2]   = coefficientsOfContinuousSpectrum('Isotropic',k_0, k_0, a_0, EE1,0,0, MU1, EE, MU, -m, (p), 0, 1);
% [B_01, B_02, DE, DH] = coefficientsOfContinuousSpectrum('Isotropic',k_0, k_0, a_0, EE1,0,0, MU1, EE, MU, m, p_0, 0, 1);
% 
%     q1 = sqrt(EE1 -  p.^2);
%     q10 = sqrt(EE1 -  p_0.^2);
%     A_1 =  1./ (k_0 * (EE1 - p.^2));
%     A_10 = 1./ (k_0 * (EE1 - p_0.^2));   
% 
%     Ez_q_0   = @(r) (B_01.* besselj(m, k_0.* r* q10));
%     Hz_q_0   = @(r) (B_02.* besselj(m, k_0.* r* q10));
%     dEz_q_0  = @(r) (B_01.* k_0.* q10 *((besselj(m, k_0.* r* q10) * m)./ (k_0.* r* q10)  - besselj(m + 1, k_0.* r* q10)));
%     dHz_q_0  = @(r) (B_02.* k_0.* q10 *((besselj(m, k_0.* r* q10) * m)./ (k_0.* r* q10)  - besselj(m + 1, k_0.* r* q10)));
%     
%     Ephi_q_plus0  = @(r)  A_10.* ((-p_0.* (m./r)).* Ez_q_0(r) + 1i * dHz_q_0(r));
%     Hrho_q_plus0  = @(r)  A_10.* ((EE1 * (m./r)).* Ez_q_0(r) - 1i * p_0 * dHz_q_0(r));
%     Hphi_q_plus0  = @(r)  A_10.* (-1i * EE1 * dEz_q_0(r) - (p_0.* (m./r)).* Hz_q_0(r));
%     Erho_q_plus0  = @(r)  A_10.* (-1i * p_0 * dEz_q_0(r) - (m./r).* Hz_q_0(r));
% 
%     
%     Ez_q   = @(r) (B_1.* besselj(-m, k_0.* r* q1));
%     Hz_q   = @(r) (B_2.* besselj(-m, k_0.* r* q1));
%     dEz_q  = @(r) (B_1.* k_0.* q1 *((besselj(-m, k_0.* r* q1) *(- m))./ (k_0.* r* q1)  - besselj(-m + 1, k_0.* r* q1)));
%     dHz_q  = @(r) (B_2.* k_0.* q1 *((besselj(-m, k_0.* r* q1) *(- m))./ (k_0.* r* q1)  - besselj(-m + 1, k_0.* r* q1)));
%     
%     Ephi_q_minus  = @(r)  A_1.* (( -p.* (-m./r)).* Ez_q(r) + 1i * dHz_q(r));
%     Hrho_q_minus  = @(r)  A_1.* ((EE1 * (-m./r)).* Ez_q(r) - 1i * p * dHz_q(r));
%     Hphi_q_minus  = @(r)  A_1.* (-1i * EE1 * dEz_q(r) - (p.* (-m./r)).* Hz_q(r));
%     Erho_q_minus  = @(r)  A_1.* (-1i * p.* dEz_q(r) - (-m./r).* Hz_q(r));
%     
%     I_into_jump = (-Ez_q(a_0).*  Hphi_q_plus0(a_0) +...
%                     Ephi_q_minus(a_0).* Hz_q_0(a_0) +...
%                     Ez_q_0(a_0).* Hphi_q_minus(a_0) -...
%                     Ephi_q_plus0(a_0).* Hz_q(a_0));
%             
%     I_into_jump = I_into_jump.* (a_0./(1i*k_0*(p + p_0)))
%     
%  
% %%% analytical formula
%     q1 = sqrt(EE1*MU1/(EE*MU) -  p.^2);
%     A_1 = 1./ (k_0 * (EE1*MU1 - EE*MU * p.^2));       
%     A = 1./ (k_0 * (1 - p.^2)); 
% %%% norm
%     G_Jm_a0 = - 1/2 * (besselj(0, k_0*q1*a_0)).^2;
%     G_Jm_0  = - 1/2 * (besselj(0, k_0*q1*0)).^2;    
%     G_Hm2_a0 = - 1/2 * (besselh(0,2, k_0*q*a_0)).^2;
%     for n=0:(m-1)
%         G_Jm_a0  =  G_Jm_a0  - (besselj(m-n, k_0*q1*a_0)).^2;
%         G_Hm2_a0 =  G_Hm2_a0 - (besselh(m-n, 2, k_0*q*a_0)).^2;
%     end
%     
%     G_Jm_1__a0 = - 1/2 * (besselj(0, k_0*q1*a_0)).^2;
%     G_Jm_1__0  = - 1/2 * (besselj(0, k_0*q1*0)).^2;
%     G_Hm2_1__a0 = - 1/2 * (besselh(0,2, k_0*q*a_0)).^2;
%     for n=0:(m-2)
%         G_Jm_1__a0  =  G_Jm_1__a0  - (besselj(m-1-n, k_0*q1*a_0)).^2;
%         G_Hm2_1__a0 =  G_Hm2_1__a0 - (besselh(m-1-n, 2, k_0*q*a_0)).^2;
%     end
%     
%     F_Jm_a0 = ((k_0*q1*a_0).^2 / 2).* ((besselj(m+1, k_0*q1*a_0)).^2  - besselj(m,   k_0*q1*a_0).*  besselj(m+2,  k_0*q1*a_0));
%     F_Hm2_a0 =((k_0*q *a_0).^2 / 2).* ((besselh(m+1,2, k_0*q*a_0)).^2 - besselh(m,2, k_0*q *a_0).* besselh(m+2,2, k_0*q *a_0));
%     
%     I_into = (-1)^m * (A_1.^2).* (-2 * p)* (EE1 * B_1_mn.* B_1_pl - B_2_mn.* B_2_pl).*...
%         (-m * (G_Jm_a0 - G_Jm_0) + m * (G_Jm_1__a0 - G_Jm_1__0) + F_Jm_a0);
%     I_into = I_into + (-1)^m * (A_1.^2) * m * 1i * (p.^2 + EE1*MU1) *...
%         ((B_1_mn.* B_2_pl + B_2_mn.* B_1_pl).* (besselj(m, k_0*q1*a_0)).^2); %%% ????????? ?????? ????????????
%     
%         
%     I_out  =  (-1)^m * (A.^2).* (-2 * p)*(Cm2_mn.* Cm2_pl - Dm2_mn.* Dm2_pl).*...
%         (-m * G_Hm2_a0 + m * G_Hm2_1__a0 + F_Hm2_a0);
%     I_out  = I_out + (-1)^m * (A.^2) * m * 1i * (p.^2 + 1) *...
%         ((Cm2_mn.* Dm2_pl + Cm2_pl.* Dm2_mn).* (besselh(m,2, k_0*q*a_0)).^2);
%     
%     N_n = (I_into - I_out);




toc
