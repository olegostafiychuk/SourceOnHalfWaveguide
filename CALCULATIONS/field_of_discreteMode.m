function component = field_of_discreteMode(componentOfField, typeOfCylinder, r, q, p, a_p_field, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, z)
  
  switch(typeOfCylinder)
      case 'Isotropic'  
          %%% ��������� ������ ��� ���� ������������ �������          
          [B_1_forward, B_2_forward, Cm2_forward,  Dm2_forward]  = coefficientsOfContinuousSpectrum(typeOfCylinder, k_0, k_0, a_0, EE1,0,0, MU1, EE, MU,  m, p, q, 0, 1);
          
          H2m  = besselh(m, 2, k_0*q*r);
          dH2m = (H2m * m)./ (k_0*q*r)  - besselh(m + 1, 2, k_0*q*r);
          A0 = 1./ (k_0 * (1 - p.^2));
          A1 = 1./ (k_0 * (EE1 - p.^2));
          q1 = sqrt(EE1*MU1 - p.^2);
          Q1 = k_0.* r * q1;
          Jm  = besselj(m, Q1);
          dJm = (Jm * m)./ Q1  - besselj(m + 1, Q1);
         
          switch(componentOfField)
              case 'Ez'
                  component = q.*(a_p_field.* Cm2_forward.* H2m).*(r>=a_0) +...
                  q1.*(a_p_field.* B_1_forward.* Jm).*(r<a_0);
              case 'Hz'
                  component = q.*(a_p_field.* Dm2_forward.* H2m).*(r>=a_0) +...
                  q1.*(a_p_field.* B_2_forward.*  Jm).*(r<a_0);
              case 'Ephi'
                  component = q.* (A0.* ((-p.* (m./r)).*...
                      (a_p_field.* Cm2_forward.*  H2m)+...
                      1i *...
                      (a_p_field.* Dm2_forward.* k_0.* q.* dH2m))).* (r>=a_0)+...
                      q1.* (A1.* ((-p.* (m./r)).*...
                      (a_p_field.* B_1_forward.*  Jm)+...
                      1i *...
                      (a_p_field.* B_2_forward.*  k_0.* q1.* dJm))).* (r<a_0);
              case 'Erho'
                  component = q.* (A0.* ((-(m./r)).*...
                      (a_p_field.* Dm2_forward.*  H2m)+...
                      -1i * p.* ...
                      (a_p_field.* Cm2_forward.*  k_0.* q.* dH2m))).* (r>=a_0)+...
                      q1.* (A1.* ((-(m./r)).*...
                      (a_p_field.* B_2_forward.*  Jm) +...
                      -1i * p.* ...
                      (a_p_field.* B_1_forward.*  k_0.* q1.* dJm))).* (r<a_0);
              case 'Hphi'
                  component = q.* (A0.* ((-p.* (m./r)).*...
                      (a_p_field.* Dm2_forward.*  H2m)+...
                      -1i *...
                      (a_p_field.* Cm2_forward.* k_0.* q.*  dH2m))).* (r>=a_0)+...
                      q1.* (A1.* ((-p.* (m./r)).*...
                      (a_p_field.* B_2_forward.* Jm) +...
                      -1i * EE1.*...
                      (a_p_field.* B_1_forward.*  k_0.* q1.* dJm))).* (r<a_0);
              case 'Hrho'
                  component = q.* (A0.* (((m./r)).*...
                      (a_p_field.* Cm2_forward.*  H2m)+...
                      -p.* 1i *...
                      (a_p_field.* Dm2_forward.* k_0.* q.* dH2m))).* (r>=a_0)+...
                      q1.* (A1.* (((EE1 * m./r)).*...
                      (a_p_field.* B_1_forward.*  Jm)+...
                      -p.* 1i *...
                      (a_p_field.* B_2_forward.*  k_0.* q1.* dJm))).* (r<a_0);
          end
      case 'Gyrotropic'
          %%% ��������� ������ ��� ���� ������������ �������          
          [B_1_forward, B_2_forward, Cm2_forward,  Dm2_forward]  =...
              coefficientsOfContinuousSpectrum(typeOfCylinder, k_0, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU,  m,   p, q, 0, 1);
          [B_1_backward,B_2_backward,Cm2_backward, Dm2_backward] =...
              coefficientsOfContinuousSpectrum(typeOfCylinder, k_0, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU,  m, (-p), q, 0, 1);
          
          H2m  = besselh(m, 2, k_0*q*r);
          dH2m = (H2m * m)./ (k_0*q*r)  - besselh(m + 1, 2, k_0*q*r);
          A0 = 1./ (k_0 * (1 - p.^2));
          
          [q1, q2, n1_back, n2_back, alp1_back, alp2_back, bet1_back, bet2_back] = term_of_gyrotropic_waveguide(EE1, GG1, HH1,-p);
          [q1, q2, n1_forw, n2_forw, alp1_forw, alp2_forw, bet1_forw, bet2_forw] = term_of_gyrotropic_waveguide(EE1, GG1, HH1, p);
          
          Q1 = (q1.* r).* k_0;
          Q2 = (q2.* r).* k_0;  
          JM1    = besselj(m+1, Q1);
          JM2    = besselj(m+1, Q2);
          Jm1    = besselj(m, Q1);
          Jm2    = besselj(m, Q2);
          
%           %%%% vary larger JM1 and Jm1
%           JM1(abs(imag(Q1))>300) = besselj(m+1, 1i*300);
%           Jm1(abs(imag(Q1))>300) = besselj(m, 1i*300);
%           %%%% vary larger JM2 and Jm2
%           JM2(abs(imag(Q2))>300) = besselj(m+1, 1i*300);
%           Jm2(abs(imag(Q2))>300) = besselj(m,   1i*300);
          
          Jm1_Q1 = Jm1./ Q1;
          Jm2_Q2 = Jm2./ Q2;
          
          a_p_field_forw = a_p_field;
          a_p_field_back = 0;
          
          switch(componentOfField)
              case 'Ez'
                  component = (a_p_field_forw.* Cm2_forward.* q.* H2m +...
                              a_p_field_back.* Cm2_backward.* q.* H2m).*(r>=a_0) +...
                HH1 * (a_p_field_forw.* (B_1_forward.*  (((1i./HH1).* n1_forw).* q1).* Jm1 +...
                                         B_2_forward.*  (((1i./HH1).* n2_forw).* q2).* Jm2)+...
                       a_p_field_back.* (B_1_backward.*  (((1i./HH1).* n1_back).* q1).* Jm1 +...
                                         B_2_backward.*  (((1i./HH1).* n2_back).* q2).* Jm2)).*(r<a_0);
              case 'Hz'
                  component = (a_p_field_forw.* Dm2_forward.*  q.* H2m +...
                               a_p_field_back.* Dm2_backward.* q.* H2m).*(r>=a_0) +...
                            (a_p_field_forw.* (B_1_forward.*  (- q1.* Jm1) +...
                                               B_2_forward.*  (- q2.* Jm2))+...
                             a_p_field_back.* (B_1_backward.*  (- q1.* Jm1) +...
                                               B_2_backward.*  (- q2.* Jm2))).*(r<a_0);
              case 'Ephi'
                  component = (A0.* ((-p.* (m./r)).*...
                      (a_p_field_forw.* Cm2_forward.*  q.* H2m+...
                      -a_p_field_back.* Cm2_backward.* q.* H2m)+...
                      1i *...
                      (a_p_field_forw.* Dm2_forward.*  k_0.* q.* q.* dH2m+...
                       a_p_field_back.* Dm2_backward.* k_0.* q.* q.* dH2m))).* (r>=a_0);
                   component = component +...
                      (a_p_field_forw.* (B_1_forward.* 1i.* (JM1 + (alp1_forw.* m).* Jm1_Q1)+...
                                         B_2_forward.* 1i.* (JM2 + (alp2_forw.* m).* Jm2_Q2))+...
                       a_p_field_back.* (B_1_backward.* 1i.* (JM1 + (alp1_back.* m).* Jm1_Q1)+...
                                         B_2_backward.* 1i.* (JM2 + (alp2_back.* m).* Jm2_Q2))).* (r<a_0);


%             [q1, q2, n1,      n2,      alp1,      alp2,      bet1,      bet2] = term_of_gyrotropic_waveguide(EE1, GG1, HH1, p);
%             [q1, q2, n1_back, n2_back, alp1_back, alp2_back, bet1_back, bet2_back] = term_of_gyrotropic_waveguide(EE1, GG1, HH1,-p);
%             
%             Q1 = (q1.* r).* k_0;
%             Q2 = (q2.* r).* k_0;  
%             
%             JM1    = besselj(m+1, Q1);
%             JM2    = besselj(m+1, Q2);
%             Jm1    = besselj(m, Q1);
%             Jm2    = besselj(m, Q2);
%             Jm1_Q1 = Jm1./ Q1;
%             Jm2_Q2 = Jm2./ Q2;
%             
%             Ez_q1   = (((1i./HH1).* n1).* q1).* Jm1;
%             Ez_q2   = (((1i./HH1).* n2).* q2).* Jm2;
%             Hz_q1   = - q1.* Jm1;
%             Hz_q2   = - q2.* Jm2;            
%             dEz_q1   = (((1i./HH1).* n1).* q1).* (m.* Jm1_Q1 - JM1).* q1.* k_0;
%             dEz_q2   = (((1i./HH1).* n2).* q2).* (m.* Jm2_Q2 - JM2).* q2.* k_0;
%             dHz_q1   = - q1.* (m.* Jm1_Q1 - JM1).* q1.* k_0;
%             dHz_q2   = - q2.* (m.* Jm2_Q2 - JM2).* q2.* k_0;
%             A_1 = 1./ (k_0 * (GG1^2 - (p.^2 - EE1).^2));
% 
%                     component = component +...
%                             A_1.*...
%                                 (a_p_field_forw.*...
%                                 (p.*(EE1-p.^2).*m./r.* (B_1_forward.* Ez_q1+ B_2_forward.* Ez_q2)+...
%                                  p.*GG1.* (B_1_forward.* dEz_q1+ B_2_forward.* dEz_q2)+...
%                                  -1i*GG1.*m./r.* (B_1_forward.* Hz_q1+ B_2_forward.* Hz_q2)+...
%                                  -1i.*(EE1-p.^2).* (B_1_forward.* dHz_q1+ B_2_forward.* dHz_q2))+...
%                                  a_p_field_back.*...
%                                 (p.*(EE1-p.^2).*m./r.* (B_1_backward.* Ez_q1+ B_2_backward.* Ez_q2)+...
%                                  p.*GG1.* (B_1_backward.* dEz_q1+ B_2_backward.* dEz_q2)+...
%                                  -1i*GG1.*m./r.* (B_1_backward.* Hz_q1+ B_2_backward.* Hz_q2)+...
%                                  -1i.*(EE1-p.^2).* (B_1_backward.* dHz_q1+ B_2_backward.* dHz_q2))).* (r<a_0);              
                                     
              case 'Erho'
                  component = (A0.* q.* ((-(m./r)).*...
                      (a_p_field_forw.* Dm2_forward.*  H2m+...
                      a_p_field_back.* Dm2_backward.* H2m)+...
                      -1i * p.* ...
                      (a_p_field_forw.* Cm2_forward.*  k_0.* q.* dH2m+...
                      -a_p_field_back.* Cm2_backward.* k_0.* q.* dH2m))).* (r>=a_0)+...
                      (a_p_field_forw.* (-B_1_forward.*  ((n1_forw.*p    + GG1)./EE1.* JM1 + (alp1_forw.* m).* Jm1_Q1)+...
                                         -B_2_forward.*  ((n2_forw.*p    + GG1)./EE1.* JM2 + (alp2_forw.* m).* Jm2_Q2))+...
                       a_p_field_back.* (-B_1_backward.* ((n1_back.*(-p) + GG1)./EE1.* JM1 + (alp1_back.* m).* Jm1_Q1)+...
                                         -B_2_backward.* ((n2_back.*(-p) + GG1)./EE1.* JM2 + (alp2_back.* m).* Jm2_Q2))).* (r<a_0);
              case 'Hphi'
                  component = q.* (A0.* ((-p.* (m./r)).*...
                      (a_p_field_forw.* Dm2_forward.*  H2m+...
                      -a_p_field_back.* Dm2_backward.* H2m)+...
                      -1i *...
                      (a_p_field_forw.* Cm2_forward.* k_0.* q.* dH2m+...
                      a_p_field_back.* Cm2_backward.* k_0.* q.* dH2m))).* (r>=a_0)+...
                      (a_p_field_forw.* (-B_1_forward.* n1_forw.* (JM1 - (bet1_forw.* m).* Jm1_Q1)+...
                                         -B_2_forward.* n2_forw.* (JM2 - (bet2_forw.* m).* Jm2_Q2))+...
                       a_p_field_back.* (-B_1_backward.* n1_back.* (JM1 - (bet1_back.* m).* Jm1_Q1)+...
                                         -B_2_backward.* n2_back.* (JM2 - (bet2_back.* m).* Jm2_Q2))).* (r<a_0);
              case 'Hrho'
                  component = (A0.* q.* (((m./r)).*...
                      (a_p_field_forw.* Cm2_forward.*  H2m+...
                      a_p_field_back.* Cm2_backward.* H2m)+...
                      -p.* 1i *...
                      (a_p_field_forw.* Dm2_forward.* k_0.* q.* dH2m+...
                      -a_p_field_back.* Dm2_backward.* k_0.* q.* dH2m))).* (r>=a_0)+...
                 -1i *(a_p_field_forw.* (B_1_forward.* (p.* JM1 - (n1_forw.* bet1_forw.* m).* Jm1_Q1)+...
                                         B_2_forward.* (p.* JM2 - (n2_forw.* bet2_forw.* m).* Jm2_Q2))+...
                       a_p_field_back.* (B_1_backward.* ((-p).* JM1 - (n1_back.* bet1_back.* m).* Jm1_Q1)+...
                                         B_2_backward.* ((-p).* JM2 - (n2_back.* bet2_back.* m).* Jm2_Q2))).* (r<a_0);
          end          
  end

component = component.* exp(-1i * k_0 * p * z);
     
     





