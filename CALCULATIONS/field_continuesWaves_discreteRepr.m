%%%%% расчёт коэффициента для волн непрерывного спектра в прямом
%%%%% направлении (в сторону среды 2)
%%%%% второе приближение в нахождении коэффициента aplus
function component = field_continuesWaves_discreteRepr(componentOfField, typeOfCylinder,r, q, p, k_0, a_0, EE1, GG1, HH1, MU1, m, z,...
          a_p_field, psi, B_1,B_2,Cm2,  Dm2)

%   p = sqrt(1-q.^2);
%   p = real(p) - 1i * abs(imag(p));
  
switch(typeOfCylinder)
    case 'Isotropic'      
        
        if(r>a_0)
            Q = k_0.* r * q;
            H2m  = besselh(m, 2, Q);
            dH2m = (H2m * m)./ Q  - besselh(m + 1, 2, Q);
            H1m  = besselh(m, 1, Q);
            dH1m = (H1m * m)./ Q  - besselh(m + 1, 1, Q);
            A0 = 1./ (k_0 * q.^2);
            switch(componentOfField)
                case 'Ez'
                    component = (a_p_field.* Cm2.*  (psi.*  H1m + H2m)).* q;
                case 'Hz'
                    component = (a_p_field.* Dm2.*  (-psi.*  H1m + H2m)).* q;
                case 'Ephi'
                    component = A0.* ((-p.* (m./r)).*...
                                (a_p_field.* Cm2.*  (psi.*  H1m + H2m))+...
                                 1i *...
                                 (a_p_field.* Dm2.*  k_0.* q.* (-psi.*  dH1m + dH2m))).* q;                    
                case 'Erho'
                     component = A0.* ((-(m./r)).*...
                           (a_p_field.* Dm2.*  (-psi.*  H1m + H2m))+...
                                -1i * p.* ...
                            (a_p_field.* Cm2.*  k_0.* q.* (psi.*  dH1m + dH2m))).* q;                        
                case 'Hphi'
                     component = A0.* ((-p.* (m./r)).*...
                           (a_p_field.* Dm2.*  (-psi.*  H1m + H2m))+...
                                -1i *...
                            (a_p_field.* Cm2.*  k_0.* q.* (psi.*  dH1m + dH2m))).* q;
                case 'Hrho'
                    component = A0.* (((m./r)).*...
                                (a_p_field.* Cm2.*  (psi.* H1m + H2m))+...
                                 -p.* 1i.*...
                                 (a_p_field.* Dm2.*  k_0.* q.* (-psi.* dH1m + dH2m))).* q; 
            end            
        else
            A1 = 1./ (k_0 * (EE1 - p.^2));
            q1 = sqrt(EE1*MU1 - p.^2);
            Q1 = k_0.* r * q1;
            Jm  = besselj(m, Q1);
            dJm = (Jm * m)./ Q1  - besselj(m + 1, Q1);
            switch(componentOfField)
                case 'Ez'
                    component = (a_p_field.* B_1.* Jm).* q1;
                case 'Hz'
                    component = (a_p_field.* B_2).* q1;
                case 'Ephi'
                    component = A1.* ((-p.* (m./r)).*...
                                (a_p_field.* B_1.*  Jm)+...
                                 1i *...
                                 (a_p_field.* B_2.*  k_0.* q1.* dJm)).* q1;
                case 'Erho'
                     component = A1.* ((-(m./r)).*...
                           (a_p_field.* B_2.* Jm)+...
                                -1i * p.* ...
                            (a_p_field.* B_1.*  k_0.* q1.* dJm)).* q1;
                    
                case 'Hphi'
                     component = A1.* ((-p.* (m./r)).*...
                           (a_p_field.* B_2.* Jm)+...
                                -1i * EE1.* ...
                            (a_p_field.* B_1.*  k_0.* q1.* dJm)).* q1;
                case 'Hrho'
                    component = A1.* (((EE1 * m./r)).*...
                                (a_p_field.* B_1.*  Jm)+...
                                -p.* 1i.*...
                                 (a_p_field.* B_2.*  k_0.* q1.* dJm)).* q1;
            end
        end

    case 'Gyrotropic'      
       
        if(r>a_0)
            Q = k_0.* r * q;
            H2m  = besselh(m, 2, Q);
            dH2m = (H2m * m)./ Q  - besselh(m + 1, 2, Q);
            H1m  = besselh(m, 1, Q);
            dH1m = (H1m * m)./ Q  - besselh(m + 1, 1, Q);
            A0 = 1./ (k_0 * q.^2);
            switch(componentOfField)
                case 'Ez'
                    component = (a_p_field.* Cm2.*  q.* (psi.*  H1m + H2m));
                case 'Hz'
                    component = (a_p_field.* Dm2.*  q.* (-psi.*  H1m + H2m));
                case 'Ephi'
                    component = A0.* q.* ((-p.* (m./r)).*...
                                (a_p_field.* Cm2.*  (psi.*  H1m + H2m))+...
                                 1i *...
                                (a_p_field.* Dm2.*  k_0.* q.* (-psi.*  dH1m + dH2m)));                    
                case 'Erho'
                     component = A0.* q.* ((-(m./r)).*...
                           (a_p_field.* Dm2.*  (-psi.*  H1m + H2m))+...
                                -1i * p.* ...
                            (a_p_field.* Cm2.*  k_0.* q.* (psi.*  dH1m + dH2m)));                        
                case 'Hphi'
                     component = A0.* q.* ((-p.* (m./r)).*...
                           (a_p_field.* Dm2.*  (-psi.*  H1m + H2m))+...
                                -1i *...
                            (a_p_field.* Cm2.*  k_0.* q.* (psi.*  dH1m + dH2m)));
                case 'Hrho'
                    component = A0.* q.* (((m./r)).*...
                                (a_p_field.* Cm2.*  (psi.* H1m + H2m))+...
                                 -p.* 1i.*...
                                 (a_p_field.* Dm2.*  k_0.* q.* (-psi.* dH1m + dH2m))); 
            end            
        else

            
            
            [q1, q2, n1,      n2,      alp1,      alp2,      bet1,      bet2] = term_of_gyrotropic_waveguide(EE1, GG1, HH1, p);
            [q1, q2, n1_back, n2_back, alp1_back, alp2_back, bet1_back, bet2_back] = term_of_gyrotropic_waveguide(EE1, GG1, HH1,-p);
            
            Q1 = (q1.* r).* k_0;
            Q2 = (q2.* r).* k_0;  
            
            JM1    = besselj(m+1, Q1);
            JM2    = besselj(m+1, Q2);
            Jm1    = besselj(m, Q1);
            Jm2    = besselj(m, Q2);
%             %%%% vary larger JM1 and Jm1
            JM1(abs(imag(Q1))>300) = besselj(m+1, 1i*300);
            Jm1(abs(imag(Q1))>300) = besselj(m, 1i*300);            
            Jm1_Q1 = Jm1./ Q1;
            Jm2_Q2 = Jm2./ Q2;
            

           
            switch(componentOfField)
                case 'Ez'
                    component =(a_p_field.* (B_1.*  (((1i./HH1).* n1).* q1).* Jm1 +...
                                             B_2.*  (((1i./HH1).* n2).* q2).* Jm2));
                case 'Hz'
                    component = (a_p_field.* (B_1.*  (- q1.* Jm1)+...
                                              B_2.*  (- q2.* Jm2)));
                             
                case 'Ephi'
                    component = 1i.*...
                                (a_p_field.* (B_1.* (JM1 + (alp1.* m).* Jm1_Q1)+...
                                              B_2.* (JM2 + (alp2.* m).* Jm2_Q2)));

                case 'Erho'
                     component = a_p_field.* (-B_1.* ((n1.*p + GG1)./EE1.* JM1 + (alp1.* m).* Jm1_Q1)+...
                                              -B_2.* ((n2.*p + GG1)./EE1.* JM2 + (alp2.* m).* Jm2_Q2));
                    
                case 'Hphi'
                     component = (a_p_field.* (-B_1.* n1.* (JM1 - (bet1.* m).* Jm1_Q1)+...
                                               -B_2.* n2.* (JM2 - (bet2.* m).* Jm2_Q2)));
                                         
                case 'Hrho'
                    component = -1i *...
                                (a_p_field.* (B_1.* (p.* JM1 - (n1.* bet1.* m).* Jm1_Q1)+...
                                              B_2.* (p.* JM2 - (n2.* bet2.* m).* Jm2_Q2)));
            end
        end
end

% component = component.* (exp( 1i * k_0 * p * z));  %%%% we consider back waves
component = component.* (exp(-1i * k_0 * p * z));  %%%% we consider back waves


