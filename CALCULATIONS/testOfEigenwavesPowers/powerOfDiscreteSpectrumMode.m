function [Pn, innerAnalytical, outerAnalytical] = powerOfDiscreteSpectrumMode(p, waveguideParameters, sourceParameters)

c = 3e10;
m = sourceParameters.m;
k_0 = sourceParameters.k_0;

EE1 = waveguideParameters.EEinner;
GG1 = waveguideParameters.GGinner;
HH1 = waveguideParameters.HHinner;
MU1 = waveguideParameters.MUinner;
EE = waveguideParameters.EEouter;
MU = waveguideParameters.MUouter;
a_0 = waveguideParameters.radius;

innerAnalytical = 0;
outerAnalytical = 0;

q = sqrt(1-p.^2);
q = real(q) - 1i * abs(imag(q));

[B_1,B_2,Cm2, Dm2] = coefficientsOfContinuousSpectrum('Gyrotropic', k_0, k_0, a_0, EE1, GG1,HH1, MU1, EE, MU,  m, p, q, 0, 1);
[q1, q2, n1, n2, alp1, alp2, bet1, bet2] = term_of_gyrotropic_waveguide(EE1, GG1, HH1, p);


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%% outer integrals %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(q == real(q) || q == 1i* imag(q))

    Q = k_0 * q * a_0;
    H2m    = besselh(m    , 2, Q);
    H2mpl1 = besselh(m + 1, 2, Q);
    H2mpl2 = besselh(m + 2, 2, Q);
                
%     Cm2 = 1;Dm2 = 1; %%%% test
%     Cm2 = 0;Dm2 = 1; %%%% test
%     Cm2 = 1;Dm2 = 0; %%%% test
    %%%%%% closed form %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    I_out_int = (a_0^2 / 8 * (-1)^(m).* (4 * H2mpl1.* H2mpl1 - 4* H2mpl2.* H2m + 8 * m / Q.^2.* H2m.* H2m));
    I_out   = - c/ (8*pi) * 2*pi *...
              (p.* (abs(Cm2).^2 + abs(Dm2).^2).* I_out_int +...
              1i* m * a_0^2./ abs(Q).^2.* (p^2.* Cm2.* conj(Dm2) - conj(Cm2).* Dm2).* abs(H2m).^2);
          
    outerAnalytical = 1;
else
    %%%%%% numerical form %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fH2m  = @(r) besselh(m, 2, k_0.* r * q);
    fdH2m = @(r) (besselh(m, 2, k_0.* r * q) * (m))./ (k_0.* r * q)  - besselh(m + 1, 2, k_0.* r * q);
        
    A = 1./ (k_0 * q.^2);  
    Ez_tq   = @(r) Cm2.* q.* fH2m(r);
    Hz_tq   = @(r) Dm2.* q.* fH2m(r);
    dEz_tq  = @(r) Cm2.* k_0.* (q).* (q).* fdH2m(r);
    dHz_tq  = @(r) Dm2.* k_0.* (q).* (q).* fdH2m(r);
        
    Ephi_q  = @(r) A.* ((-p.* (m./r)).* Ez_tq(r) + 1i * dHz_tq(r));
    Hphi_q  = @(r) A.* (-1i * dEz_tq(r) - (p.* (m./r)).* Hz_tq(r));
    Hrho_q  = @(r) A.* (((m./r)).* Ez_tq(r) - 1i * p * dHz_tq(r));
    Erho_q  = @(r) A.* (-1i * p * dEz_tq(r) - (m./r).* Hz_tq(r));
    
    I_out = c/ (8*pi) * 2*pi *...
        integral(@(r) (Erho_q(r).* conj(Hphi_q(r)) - Ephi_q(r).* conj(Hrho_q(r))).* r, a_0, Inf);
        
end
%%%%%%%%%%%%%%%%%%% end outer integrals %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% inner integrals %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if((q1 == real(q1) || q1 == 1i* imag(q1)) && (q2 == real(q2) || q2 == 1i* imag(q2)))
%     B_1 = 1;B_2 = 1; %%%% test
%     B_1 = 1;B_2 = 0; %%%% test
%     B_1 = 0;B_2 = 1; %%%% test
%      q1 = 10; q2 = -1i * 10; %%%%% test
%      q2 = 10; q1 = -1i * 10; %%%%% test
    %%%%%% closed form %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    JM1    = besselj(m+1, (q1.* a_0).* k_0);
    JM2    = besselj(m+1, (q2.* a_0).* k_0);
    JMM1    = besselj(m+2, (q1.* a_0).* k_0);
    JMM2    = besselj(m+2, (q2.* a_0).* k_0);
    Jm1    = besselj(m, (q1.* a_0).* k_0);
    Jm2    = besselj(m, (q2.* a_0).* k_0);
    Q1 = q1.* a_0.* k_0;
    Q2 = q2.* a_0.* k_0;
    
    gamma_1 = q1./ conj(q1); gamma_2 = q2./ conj(q2);
    I1_1 = (a_0.^2 / 2).* (JM1.^2  - Jm1.*  JMM1).* gamma_1.^(m+1);
    I1_2 = (a_0.^2 / 2).* (JM2.^2  - Jm2.*  JMM2).* gamma_2.^(m+1);
    I5_12 =  (a_0.^2./ (Q1^2 - Q2^2)).*...
                 (Q1.* JMM1.* JM2 - Q2.* JMM2.* JM1).* gamma_2.^(m+1);
    I5_21 =  (a_0.^2./ (Q1^2 - Q2^2)).*...
                 (Q1.* JMM1.* JM2 - Q2.* JMM2.* JM1).* gamma_1.^(m+1);
             
    I_in = c/ (8*pi) * 2*pi.*... 
                    (abs(B_1)^2.* (I1_1.* (p - conj(n1).* (1 + alp1)) -...
                                   m./ (k_0* abs(q1))^2 * alp1.* conj(bet1).* conj(n1).* Jm1.^2.* gamma_1.^m) +...
                     abs(B_2)^2.* (I1_2.* (p - conj(n2).* (1 + alp2)) -...
                                   m./ (k_0* abs(q2))^2 * alp2.* conj(bet2).* conj(n2).* Jm2.^2.* gamma_2.^m) +...
                     B_1.* conj(B_2).* (I5_12.* (p - conj(n2).* (1 + alp1)) +...
                     - m.* a_0^2/ (Q1.* conj(Q2)).* alp1.* conj(bet2).* conj(n2).* Jm1.* Jm2.* gamma_2.^m) +...
                     B_2.* conj(B_1).* (I5_21.* (p - conj(n1).* (1 + alp2)) +...
                     - m.* a_0^2/ (Q2.* conj(Q1)).* alp2.* conj(bet1).* conj(n1).* Jm1.* Jm2.* gamma_1.^m));
                 
     innerAnalytical = 1;    
else
    %%%%%% numerical form %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
          JM1    = @(r) besselj(m+1, (q1.* r).* k_0);
          JM2    = @(r) besselj(m+1, (q2.* r).* k_0);
          Jm1    = @(r) besselj(m, (q1.* r).* k_0);
          Jm2    = @(r) besselj(m, (q2.* r).* k_0);
          Jm1_Q1 = @(r) Jm1(r)./ ((q1.* r).* k_0);
          Jm2_Q2 = @(r) Jm2(r)./ ((q2.* r).* k_0);
          
%           Ez_q    = @(r) B_1.*  (((1i./HH1).* n1).* q1).* Jm1(r) +...
%                          B_2.*  (((1i./HH1).* n2).* q2).* Jm2(r);
%           Hz_q    = @(r) - B_1.*  q1.* Jm1(r) +...
%                          - B_2.*  q2.* Jm2(r);
          Ephi_q  = @(r)  B_1 * 1i * (JM1(r) + (alp1.* m).* Jm1_Q1(r))+...
                          B_2 * 1i * (JM2(r) + (alp2.* m).* Jm2_Q2(r));
          Hrho_q  = @(r) - 1i * B_1 * (p*JM1(r) - n1.* (bet1 * m).* Jm1_Q1(r))+...
                         - 1i * B_2 * (p*JM2(r) - n2.* (bet2 * m).* Jm2_Q2(r));
          Hphi_q  = @(r) - B_1 * n1.* (JM1(r) - (bet1 * m).* Jm1_Q1(r))+...
                         - B_2 * n2.* (JM2(r) - (bet2 * m).* Jm2_Q2(r));
          Erho_q  = @(r) - B_1 * ((n1*(p) + GG1)/EE1 * JM1(r) + (alp1.* m).* Jm1_Q1(r))+...
                         - B_2 * ((n2*(p) + GG1)/EE1 * JM2(r) + (alp2.* m).* Jm2_Q2(r));

    I_in = c/ (8*pi) * 2*pi *...
         integral(@(r) (Erho_q(r).* conj(Hphi_q(r)) - Ephi_q(r).* conj(Hrho_q(r))).* r, 0, a_0,'RelTol',1e-8,'AbsTol',1e-13);
end
%%%%%%%%%%%%%%%%%%% end inner integrals %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   
 Pn = real(I_in + I_out);


