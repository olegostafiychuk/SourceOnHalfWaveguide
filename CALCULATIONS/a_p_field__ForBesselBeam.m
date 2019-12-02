%%%%% расчёт коэффициента для волн непрерывного спектра в праямом
%%%%% направлении (в сторону среды 2)
%%%%% второе приближение в нахождении коэффициента aplus
function aplus2 = a_p_field__ForBesselBeam(typeOfCylinder, q, q_0, p, k_0, k, a_0, EE1, GG1, HH1, MU1, EE2, MU2, EE, MU, m, z, psi, AE_0, AH_0)

c = 3e10;   % velocity of light

p_0 = sqrt(1-q_0.^2);
p_0 = real(p_0) - 1i * abs(imag(p_0));


% m0 = abs(m);
m0 = (m);

switch(typeOfCylinder)
    case 'Isotropic'
        m0 = (m);
        m = -m;
        p = -p;        
        [B_1, B_2, Cm2, Dm2]  = coefficientsOfContinuousSpectrum(typeOfCylinder, k_0, k, a_0, EE1, GG1, HH1, MU1, EE, MU, m, p, q, psi,  1);
        %%%%%% norm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        N_p = (-1)^m * 8 * (p./ q / k_0).* (Cm2.^2 + Dm2.^2).* psi;
        
        % Dm2=1;
        %%% нормированные коэффициенты возбуждения волн
        % psi = psi1_q(k_0, k, a_0, EE1, MU1, EE, MU, m, p);
        % [B_1,B_2,Cm2, Dm2] = coefficientsOfContinuousSpectrum(k_0, k, a_0, EE1, MU1, EE, MU, m, p, q, psi, Dm2);
        
        
        
        % r = a_0;
        %%% структура поля моды
        H2m  = besselh(m, 2, k_0.* a_0 * q);
        dH2m = (besselh(m, 2, k_0.* a_0 * q) * m)./ (k_0.* a_0 * q)  - besselh(m + 1, 2, k_0.* a_0 * q);
        H1m  = besselh(m, 1, k_0.* a_0 * q);
        dH1m = (besselh(m, 1, k_0.* a_0 * q) * m)./ (k_0.* a_0 * q)  - besselh(m + 1, 1, k_0.* a_0 * q);
        
        A = 1./ (k_0 * (1 - p.^2));  
        Ez_q  = q.* Cm2.* (psi.* H1m + H2m);
        Hz_q  = q.* Dm2.* (-psi.* H1m + H2m);
        Ephi_q  = A.* ((-p.* (m./a_0)).* Ez_q + 1i * k_0 * q.* q.*...
            Dm2.* (-psi.* dH1m + dH2m));
        Hphi_q  = A.* ((-p.* (m./a_0)).* Hz_q - 1i * k_0 * q.* q.*...
            Cm2.* (psi.* dH1m + dH2m));
        
        %%%%% type Of Bessel Beam is H-polarized beam
        A0 = 1./ (k_0 * (1 - p_0.^2));
        %   Hz_inc    = AE_0.* (besselj(m0, k_0.* a_0 * q_0));                       
        %   Ez_inc    = 0;
        %   Hphi_inc  = AE_0.* ((k_0 * q_0.^2)^(-1)).* (- (p_0.* (m0./a_0)).*  (besselj(m0, k_0.* a_0 * q_0)));
        %   Ephi_inc  = AE_0.* ((k_0 * q_0.^2)^(-1)).* (1i * MU * k_0 * (q_0).*...
        %       (besselj(m0-1, k_0.* a_0 * q_0) - besselj(m0+1, k_0.* a_0 * q_0))/2);
        %   [B_01,B_02,DE,DH] = singleCylinderCoefficients_BesselBeam_forTest(k_0 * c, a_0, EE1, MU1, EE, MU,  c, 1, m0, p_0, q_0, 0, Ez_inc, Hz_inc, Ephi_inc, Hphi_inc);
        
        %%%%% type Of Bessel Beam is H-polarized beam
        %     Ez_q_0   = 0;
        %     Hz_q_0   = AE_0.* besselj(m0,  k_0.* a_0* q_0);
        %     dEz_q_0  = 0;
        %     dHz_q_0  = AE_0.* k_0.* (q_0).*((besselj(m0,   k_0.* a_0* q_0) * (m0))./ (k_0.* a_0* q_0)  - besselj(m0 + 1,   k_0.* a_0* q_0));
        %     
        %     Ephi_q_0  = A0.* ((-p_0.* (m0./a_0)).* Ez_q_0 + 1i * dHz_q_0);
        %     Hphi_q_0  = A0.* (-1i * dEz_q_0 - (p_0.* (m0./a_0)).* Hz_q_0);
        
        Ez_q_0   = @(r) AE_0.* besselj(m0,  k_0.* r* q_0);
        Hz_q_0   = @(r) AH_0.* besselj(m0,  k_0.* r* q_0);
        dEz_q_0  = @(r) AE_0.* k_0.* (q_0) *((besselj(m0, k_0.* r* q_0) * (m0))./ (k_0.* r* q_0)  - besselj(m0 + 1,   k_0.* r* q_0));
        dHz_q_0  = @(r) AH_0.* k_0.* (q_0) *((besselj(m0, k_0.* r* q_0) * (m0))./ (k_0.* r* q_0)  - besselj(m0 + 1,   k_0.* r* q_0));
        
        Ephi_q_plus0  = @(r)  A0.* ((-p_0.* (m0./r)).* Ez_q_0(r) + 1i * dHz_q_0(r));
        Hphi_q_plus0  = @(r)  A0.* (-1i * dEz_q_0(r) - (p_0.* (m0./r)).* Hz_q_0(r));
        %     Hrho_q_plus0  = @(r)  A0.* (((m0./r)).* Ez_q_0(r) - 1i * p_0 * dHz_q_0(r));
        %     Erho_q_plus0  = @(r)  A0.* (-1i * p_0 * dEz_q_0(r) - (m0./r).* Hz_q_0(r));
        
        %%%% из теоремы Гаусса
        S_plus = (-Ez_q.* Hphi_q_plus0(a_0) + Ephi_q.* Hz_q_0(a_0) +...
            Ez_q_0(a_0).* Hphi_q - Ephi_q_plus0(a_0).* Hz_q).* (a_0 * exp(-1i*(p+p_0)*k_0*z)./(1i*k_0*(p+p_0)));
        
        
        
        
        
        %      H2m = @(r)  besselh(m, 2, k_0.* r * q);
        %     dH2m = @(r) k_0.* q * ((besselh(m, 2, k_0.* r * q) * m)./ (k_0.* r * q)  - besselh(m + 1, 2, k_0.* r * q));
        %     H1m  = @(r) besselh(m, 1, k_0.* r * q);
        %     dH1m = @(r) k_0.* q * ((besselh(m, 1, k_0.* r * q) * m)./ (k_0.* r * q)  - besselh(m + 1, 1, k_0.* r * q));
        %         
        %        A = 1./ (k_0 * (1 - p.^2));  
        %       Ez_q  = @(r)  Cm2.* (psi.* H1m(r) + H2m(r));
        %       Hz_q  = @(r)  Dm2.* (-psi.* H1m(r) + H2m(r));
        %       dEz_q = @(r)  Cm2.* (psi.* dH1m(r) + dH2m(r));
        %       dHz_q = @(r)  Dm2.* (-psi.* dH1m(r) + dH2m(r));
        %     
        %     Ephi_q_out_minus  = @(r)  A.* ((-p.* (m./r)).* Ez_q(r) + 1i * dHz_q(r));
        %     Hrho_q_out_minus  = @(r)  A.* (((m./r)).* Ez_q(r) - 1i * p * dHz_q(r));
        %     Hphi_q_out_minus  = @(r)  A.* (-1i * dEz_q(r) - (p.* (m./r)).* Hz_q(r));
        %     Erho_q_out_minus  = @(r)  A.* (-1i * p * dEz_q(r) - (m./r).* Hz_q(r));
        %    
        %     
        %     %%%%% type Of Bessel Beam is H-polarized beam
        %     Ez_q0   = @(r) 0;
        %     Hz_q0   = @(r) AE_0.* besselj(m0,  k_0.* r* q_0);
        %     dEz_q0  = @(r) 0;
        %     dHz_q0  = @(r) AE_0.* k_0.* (q_0) *((besselj(m0,   k_0.* r* q_0) * (m0))./ (k_0.* r* q_0)  - besselj(m0 + 1,   k_0.* r* q_0));
        %     
        %     Ephi_q_out_plus0  = @(r)  A0.* ((-p_0.* (m0./r)).* Ez_q0(r) + 1i * dHz_q0(r));
        %     Hphi_q_out_plus0  = @(r)  A0.* (-1i * dEz_q0(r) - (p_0.* (m0./r)).* Hz_q0(r));
        %     Hrho_q_out_plus0  = @(r)  A0.* (((m0./r)).* Ez_q0(r) - 1i * p_0 * dHz_q0(r));
        %     Erho_q_out_plus0  = @(r)  A0.* (-1i * p_0 * dEz_q0(r) - (m0./r).* Hz_q0(r));
        %     
        %     
        %     S_plus = integral(@(r) (Erho_q_out_plus0(r).* Hphi_q_out_minus(r) - Ephi_q_out_plus0(r).* Hrho_q_out_minus(r)-...
        %         Erho_q_out_minus(r).* Hphi_q_out_plus0(r) + Ephi_q_out_minus(r).* Hrho_q_out_plus0(r)).* r, a_0, 100 * a_0)
        %     S_plus =-(-Ez_q(a_0).* Hphi_q_out_plus0(a_0) + Ephi_q_out_minus(a_0).* Hz_q0(a_0) +...
        %                Ez_q0(a_0).* Hphi_q_out_minus(a_0) - Ephi_q_out_plus0(a_0).*  Hz_q(a_0)).* (a_0./(1i*k_0*(p + p_0)))
        
        
        
        
        
        
        
        
        %     %%% структура поля моды
        q1 = sqrt(EE1*MU1 - p.^2);
        A_1 = 1./ (k_0 * (EE1*MU1 - p.^2));
        Jm  = @(r) besselj(m, k_0.* r * q1);
        dJm = @(r) (Jm(r) * m)./ (k_0.* r * q1)  - besselj(m + 1, k_0.* r * q1);
        
        Ez_q  = @(r) q1.* B_1.* Jm(r);
        Hz_q  = @(r) q1.* B_2.* Jm(r);
        Ephi_q  = @(r) A_1.* ((-p.* (m./r)).* Ez_q(r) + 1i * k_0 * q1.* q1.*...
            B_2.* dJm(r));
        Hphi_q  = @(r) A_1.* ((-p.* (m./r)).* Hz_q(r) - 1i * k_0 * EE1 * q1.* q1.*...
            B_1.* dJm(r));    
        Hrho_q  = @(r)  A_1.* ((EE1 * (m./r)).* Ez_q(r) - 1i * k_0 * p.* q1.* q1.*...
            B_2.* dJm(r));
        Erho_q  = @(r)  A_1.* (-1i * k_0 * p.* q1.* q1.*...
            B_1.* dJm(r) - (m./r).* Hz_q(r));
        
        % %%% структура тестового поля
        q10 = sqrt(EE1*MU1 -  p_0.^2);
        A_0 = 1./ (k_0 * (1 - p_0.^2));
        Jm_0  = @(r) besselj(m0, k_0.* r * q_0);
        dJm_0 = @(r) (Jm_0(r) * m0)./ (k_0.* r * q_0)  - besselj(m0 + 1, k_0.* r * q_0);
        
        %     Ez_q_0  = @(r) 0;
        %     Hz_q_0  = @(r) AE_0.* Jm_0(r);
        %     Ephi_q_0 =  @(r) A_0.* ( + 1i * k_0 * q_0.*...
        %         AE_0.* dJm_0(r));
        %     Hphi_q_0 = @(r)  A_0.* ((-p_0.* (m0./r)).* Hz_q_0(r));
        %     Hrho_q_0  = @(r)  A_0.* ( - 1i * k_0 * p_0.* q_0.*...
        %         AE_0.* dJm_0(r));
        %     Erho_q_0  = @(r)  A_0.* (- (m0./r).* Hz_q_0(r));    
        
        %     S_minus = integral(@(r) (Erho_q_0(r).* Hphi_q(r) - Ephi_q_0(r).* Hrho_q(r)-...
        %         Erho_q(r).* Hphi_q_0(r) + Ephi_q(r).* Hrho_q_0(r)).* r, 0, a_0)
        
        
        %     Jm0  = @(r) besselj(m0,  k_0.* r* q_0);
        %     dJm0 = @(r) k_0.* q_0 *((besselj(m0,   k_0.* r* q_0) * (m0))./ (k_0.* r* q_0)  - besselj(m0 + 1,   k_0.* r* q_0));
        %     Jm   = @(r) besselj(m0,  k_0.* r* q1);
        %     dJm  = @(r) k_0.* q1 *((besselj(m0,   k_0.* r* q1) * (m0))./ (k_0.* r* q1)  - besselj(m0 + 1,   k_0.* r* q1));
        
        %     S_minus = (-1)^m0 * integral(@(r) A_1.* A0.* AE_0 * ((1i * (m0./r).* (EE1 - p.* p_0)).* B_1.* (Jm0(r).* dJm(r) + dJm0(r).* Jm(r)) +...
        %         ((p_0 - p) * B_2).* ((m0^2./ r.^2).* Jm0(r).* Jm(r) + dJm0(r).* dJm(r))).* r, 0, a_0)
        
        
        %     S_minus = integral(@(r) (Erho_q_plus0(r).* Hphi_q_minus(r) - Ephi_q_plus0(r).* Hrho_q_minus(r)-...
        %         Erho_q_minus(r).* Hphi_q_plus0(r) + Ephi_q_minus(r).* Hrho_q_plus0(r)).* r, 0, a_0);
        
        I02 = 1/2*((a_0)./ ((k_0*q_0).^2 - (k_0*q1).^2)).* (k_0.* q_0.* besselj(m0  , k_0.* a_0* q_0).* besselj(m0-1, k_0.* a_0* q1)  - k_0.* q1.* besselj(m0-1, k_0.* a_0* q_0).* besselj(m0  , k_0.* a_0* q1))-...
            1/2*((a_0)./ ((k_0*q_0).^2 - (k_0*q1).^2)).* (k_0.* q_0.* besselj(m0+2, k_0.* a_0* q_0).* besselj(m0+1, k_0.* a_0* q1)  - k_0.* q1.* besselj(m0+1, k_0.* a_0* q_0).* besselj(m0+2, k_0.* a_0* q1));
        I03 = ((a_0)./ ((k_0*q_0).^2 - (k_0*q1).^2)).* (k_0.* q_0.* besselj(m0+2, k_0.* a_0* q_0).* besselj(m0+1, k_0.* a_0* q1)  - k_0.* q1.* besselj(m0+1, k_0.* a_0* q_0).* besselj(m0+2, k_0.* a_0* q1));
        I04 =   besselj(m0, k_0.* a_0* q_0).* besselj(m0, k_0.* a_0* q1);
        S_minus = (-1)^m0 * (A_1.* A_0.* AH_0.* ((p_0 - p).* B_2).* I02.* (k_0.*q_0).* (k_0.*q1)+...
            A_1.* A_0.* AH_0.* ((p_0 - p).* B_2).* I03.* (k_0.*q_0).* (k_0.*q1)+...
            A_1.* A_0.* AH_0.* (1i * (m0).* (EE1 - p.* p_0)).* B_1.* I04);
        
        % S_minus = -(-1)^m0 * A_1.* A_0.* AE_0.* (((p_0*EE1 - p).* B_1).* I02.* (k_0.*q_0).* (k_0.*q1)+...
        %                                          ((p_0*EE1 - p).* B_1).* I03.* (k_0.*q_0).* (k_0.*q1)+...
        %                                          (1i * (m0).* (p.* p_0 - 1)).* B_2.* I04+...
        %                                          (1i * (m0).* (p.* p_0 - EE1)).* B_1.* I04);
        
        S_minus = S_minus-(-1)^m0 * A_1.* A_0.* AE_0.* ((((p_0*EE1 - p).* B_1)).* (I03.* (k_0.*q_0).* (k_0.*q1) + m0 * I04)+...
                                          ((1i * (m0).* (p.* p_0 - 1)).* B_2).* I04);
                                      
        S_minus = S_minus.* q1;
                 
             
                 
                 
   
%     q1 = sqrt(EE1*MU1 -  p.^2);
%     A_1 = 1./ (k_0 * (EE1*MU1 - EE*MU * p.^2));
%     Ez_q   = @(r) B_1.* besselj(m, k_0.* r* q1);
%     Hz_q   = @(r) B_2.* besselj(m, k_0.* r* q1);
%     dEz_q  = @(r) B_1.* k_0.* q1.*((besselj(m, k_0.* r* q1) * m)./ (k_0.* r* q1)  - besselj(m + 1, k_0.* r* q1));
%     dHz_q  = @(r) B_2.* k_0.* q1.*((besselj(m, k_0.* r* q1) * m)./ (k_0.* r* q1)  - besselj(m + 1, k_0.* r* q1));
%     Ephi_q_minus  = @(r)  A_1.* (( -p.* (m./r)).* Ez_q(r) + 1i * dHz_q(r));
%     Hrho_q_minus  = @(r)  A_1.* ((EE1 * (m./r)).* Ez_q(r) - 1i * p.* dHz_q(r));
%     Hphi_q_minus  = @(r)  A_1.* (-1i * EE1 * dEz_q(r) - (p.* (m./r)).* Hz_q(r));
%     Erho_q_minus  = @(r)  A_1.* (-1i * p.* dEz_q(r) - (m./r).* Hz_q(r));            
%                 
%     Ez_q0   = @(r) AE_0.* besselj(m0,  k_0.* r* q_0);
%     Hz_q0   = @(r) AH_0.* besselj(m0,  k_0.* r* q_0);
%     dEz_q0  = @(r) AE_0.* k_0.* (q_0).*((besselj(m0,   k_0.* r* q_0) * (m0))./ (k_0.* r* q_0)  - besselj(m0 + 1,   k_0.* r* q_0));
%     dHz_q0  = @(r) AH_0.* k_0.* (q_0).*((besselj(m0,   k_0.* r* q_0) * (m0))./ (k_0.* r* q_0)  - besselj(m0 + 1,   k_0.* r* q_0));
%     Ephi_q_plus0  = @(r)  A0.* ((-p_0.* (m0./r)).* Ez_q0(r) + 1i * dHz_q0(r));
%     Hphi_q_plus0  = @(r)  A0.* (-1i * dEz_q0(r) - (p_0.* (m0./r)).* Hz_q0(r));
%     Hrho_q_plus0  = @(r)  A0.* (((m0./r)).* Ez_q0(r) - 1i * p_0 * dHz_q0(r));
%     Erho_q_plus0  = @(r)  A0.* (-1i * p_0 * dEz_q0(r) - (m0./r).* Hz_q0(r));
% 
%     S_minus = quadv(@(r) (Erho_q_plus0(r).* Hphi_q_minus(r) - Ephi_q_plus0(r).* Hrho_q_minus(r)-...
%         Erho_q_minus(r).* Hphi_q_plus0(r) + Ephi_q_minus(r).* Hrho_q_plus0(r)).* r, 0, a_0);    
                 
                 
%     S_minus = (-Ez_q(a_0).* Hphi_q_0(a_0) + Ephi_q(a_0).* Hz_q_0(a_0) +...
%         Ez_q_0(a_0).* Hphi_q(a_0) - Ephi_q_0(a_0).* Hz_q(a_0)).* (a_0 * exp(-1i*(p+p_0)*k_0*z)./(1i*k_0*(p+p_0)))
   
%   %%% структура поля моды
%     q1 = sqrt(EE1*MU1/(EE*MU) -  p.^2);
%     A_1 = 1./ (k_0 * (EE1*MU1 - EE*MU * p.^2));
%     Jm  = besselj(m, k_0.* r * q1);
%     dJm = (Jm * m)./ (k_0.* r * q1)  - besselj(m + 1, k_0.* r * q1);
%    
%     Ez_q  = B_1.* Jm;
%     Hz_q  = B_2.* Jm;
%     Ephi_q  = A_1.* ((-p.* (m./r)).* Ez_q + 1i * k_0 * q1.*...
%         B_2.* dJm);
%     Hphi_q  = A_1.* ((-p.* (m./r)).* Hz_q - 1i * k_0 * EE1 * q1.*...
%         B_1.* dJm);
%     
% %%%%% type Of Bessel Beam is H-polarized beam
%     q10 = sqrt(EE1*MU1/(EE*MU) -  p_0.^2);
%     A_10 = 1./ (k_0 * (EE1*MU1 - EE*MU * p_0.^2));
% 
%     Ez_q_0   = 0;
%     Hz_q_0   = AE_0.* besselj(m0,  k_0.* r* q_0);
%     dEz_q_0  = 0;
%     dHz_q_0  = AE_0.* k_0.* (q_0).*((besselj(m0,   k_0.* r* q_0) * (m0))./ (k_0.* r* q_0)  - besselj(m0 + 1,   k_0.* r* q_0));
%     
%     Ephi_q_0  = A0.* ((-p_0.* (m0./r)).* Ez_q_0 + 1i * dHz_q_0);
%     Hphi_q_0  = A0.* (-1i * dEz_q_0 - (p_0.* (m0./r)).* Hz_q_0);
% 
% 
%     S_minus = (-Ez_q.* Hphi_q_0 + Ephi_q.* Hz_q_0 +...
%         Ez_q_0.* Hphi_q - Ephi_q_0.* Hz_q).* (a_0 * exp(-1i*(p+p_0)*k_0*z)./(1i*k_0*(p+p_0)));
 
    
         aplus2 = (S_plus + S_minus)./ N_p; %%% ????? ?? ????? 1i?
    
         
         
    case 'Gyrotropic'
%         
%          p= -9.747840837847118*1i;                               
%          psi = psi1_q__gyrotropic(k_0, k, a_0, EE1, -GG1, HH1, -m, (-p));
        
        [q1, q2, n1, n2, alp1, alp2, bet1, bet2] = term_of_gyrotropic_waveguide(EE1, -GG1, HH1, -p);
        [B_1, B_2, Cm2, Dm2]  = coefficientsOfContinuousSpectrum(typeOfCylinder, k_0, k, a_0, EE1, -GG1, HH1, MU1, EE, MU, -m, -p, q, psi,  1);
        %%%%%% norm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        N_p = -(-1)^m * 8 * (p./ q / k_0).* (Cm2.^2 + Dm2.^2).* psi;
        
        % r = a_0;
        %%% структура поля моды
        H2m  = besselh(m, 2, k_0.* a_0 * q);
        dH2m = (besselh(m, 2, k_0.* a_0 * q) * (m))./ (k_0.* a_0 * q)  - besselh(m + 1, 2, k_0.* a_0 * q);
        H1m  = besselh(m, 1, k_0.* a_0 * q);
        dH1m = (besselh(m, 1, k_0.* a_0 * q) * (m))./ (k_0.* a_0 * q)  - besselh(m + 1, 1, k_0.* a_0 * q);
        
        A = 1./ (k_0 * (1 - p.^2));  
        Ez_q  = Cm2.* q.* (psi.* H1m + H2m);
        Hz_q  = Dm2.* q.* (-psi.* H1m + H2m);
        Ephi_q  = A.* ((-p.* (m./a_0)).* Ez_q + 1i * k_0 * q.*...
            Dm2.* q.* (-psi.* dH1m + dH2m));
        Hphi_q  = A.* ((-p.* (m./a_0)).* Hz_q - 1i * k_0 * q.*...
            Cm2.* q.* (psi.* dH1m + dH2m));
        
        %%%%% type Of Bessel Beam is H-polarized beam
        A0 = 1./ (k_0 * (1 - p_0.^2));
        
        Ez_q_0   = @(r) AE_0.* besselj(m0,  k_0.* r* q_0);
        Hz_q_0   = @(r) AH_0.* besselj(m0,  k_0.* r* q_0);
        dEz_q_0  = @(r) AE_0.* k_0.* (q_0) *((besselj(m0, k_0.* r* q_0) * (m0))./ (k_0.* r* q_0)  - besselj(m0 + 1,   k_0.* r* q_0));
        dHz_q_0  = @(r) AH_0.* k_0.* (q_0) *((besselj(m0, k_0.* r* q_0) * (m0))./ (k_0.* r* q_0)  - besselj(m0 + 1,   k_0.* r* q_0));
        
        Ephi_q_plus0  = @(r)  A0.* ((-p_0.* (m0./r)).* Ez_q_0(r) + 1i * dHz_q_0(r));
        Hphi_q_plus0  = @(r)  A0.* (-1i * dEz_q_0(r) - (p_0.* (m0./r)).* Hz_q_0(r));
        %     Hrho_q_plus0  = @(r)  A0.* (((m0./r)).* Ez_q_0(r) - 1i * p_0 * dHz_q_0(r));
        %     Erho_q_plus0  = @(r)  A0.* (-1i * p_0 * dEz_q_0(r) - (m0./r).* Hz_q_0(r));
        
        %%%% из теоремы Гаусса
        S_plus = (-1)^m * (-Ez_q.* Hphi_q_plus0(a_0) + Ephi_q.* Hz_q_0(a_0) +...
                           Ez_q_0(a_0).* Hphi_q - Ephi_q_plus0(a_0).* Hz_q).*...
                          (a_0 * exp(-1i*(-p+p_0)*k_0*z)./(1i*k_0*(-p+p_0)));
    
          
          
%       %     %%% структура поля моды
          JM0    = besselj(m+1, (q_0.* a_0).* k_0);                            
          JM1    = besselj(m+1, (q1.* a_0).* k_0);
          JM2    = besselj(m+1, (q2.* a_0).* k_0);
          JMM0   = besselj(m+2, (q_0.* a_0).* k_0);
          JMM1   = besselj(m+2, (q1.* a_0).* k_0);
          JMM2   = besselj(m+2, (q2.* a_0).* k_0);
          Jm0    = besselj(m, (q_0.* a_0).* k_0);
          Jm1    = besselj(m, (q1.* a_0).* k_0);
          Jm2    = besselj(m, (q2.* a_0).* k_0);
          %%%% vary larger JM1 and Jm1
          JMM1(abs(imag(q1.* a_0.* k_0))>300) = besselj(m+2, 1i*300);
          JM1(abs(imag(q1.* a_0.* k_0))>300) = besselj(m+1, 1i*300);
          Jm1(abs(imag(q1.* a_0.* k_0))>300) = besselj(m,   1i*300);
          
          M10 =(a_0./ ((k_0*q1).^2 - (k_0*q_0).^2)).*...
              (k_0*q1.* JMM1.* JM0 - k_0*q_0.* JMM0.* JM1);
          M20 =(a_0./ ((k_0*q2).^2 - (k_0*q_0).^2)).*...
              (k_0*q2.* JMM2.* JM0 - k_0*q_0.* JMM0.* JM2);
          
          A_0 = 1./ (k_0 * (1 - p_0.^2));
          A_1 = 1./ (k_0 * (GG1^2 - (p.^2 - EE1).^2));
          
          I_into_1 = (p.* (EE1 - p.^2) - p_0* (GG1^2 - EE1.* (EE1 - p.^2)))*...
              (1i/HH1)* AE_0.*...
              ((B_1.* n1.* q1).* (m.* Jm1.* Jm0 + k_0^2 * q1.* q_0.* M10)+...
               (B_2.* n2.* q2).* (m.* Jm2.* Jm0 + k_0^2 * q2.* q_0.* M20));
          
          I_into_2 = ((EE1 - p.^2).* (-p - p_0))*...
              AH_0.*...
              ((-B_1.* q1).* (m.* Jm1.* Jm0 + k_0^2 * q1.* q_0.* M10)+...
               (-B_2.* q2).* (m.* Jm2.* Jm0 + k_0^2 * q2.* q_0.* M20));
           
          I_into_3 = m * (p.* p_0 + 1).* p.* GG1.*...
                  ((1i./HH1).* (B_1.* n1.* q1.* Jm1 + B_2.* n2.* q2.* Jm2)).*...
                               (AE_0.* Jm0);
           
          I_into_4 = GG1 * (-p - p_0).*...
              m.* (-B_1.* q1.* Jm1 - B_2.* q2.* Jm2).* (AH_0.* Jm0);
          
          I_into_5 = 1i * (-p.* p_0 - 1).* (EE1 - p.^2).*...
              m.* (-B_1.* q1.* Jm1 - B_2.* q2.* Jm2).* (AE_0.* Jm0);
          
          I_into_6 = 1i * (-p.* p_0.* (EE1 - p.^2) + (GG1^2 - EE1.* (EE1 - p.^2))).*...
              m.* ((1i./HH1).* (B_1.* n1.* q1.* Jm1 + B_2.* n2.* q2.* Jm2)).* (AH_0.* Jm0);

          I_into_7 = 1i * GG1 * (-p.*p_0 - 1).*...
               AE_0.*...
              ((-B_1.* q1).* (m.* Jm1.* Jm0 + k_0^2 * q1.* q_0.* M10)+...
               (-B_2.* q2).* (m.* Jm2.* Jm0 + k_0^2 * q2.* q_0.* M20));
           
           I_into_8 = 1i * GG1 * p.* (-p_0 - p).*...
              (1i/HH1)* AH_0.*...
              ((B_1.* n1.* q1).* (m.* Jm1.* Jm0 + k_0^2 * q1.* q_0.* M10)+...
               (B_2.* n2.* q2).* (m.* Jm2.* Jm0 + k_0^2 * q2.* q_0.* M20));
          
          S_minus = (-1)^m * (A_1.* A_0).* (I_into_1 + I_into_2 + I_into_3 +...
                                            I_into_4 + I_into_5 + I_into_6 +...
                                            I_into_7 + I_into_8);
                                        
                                         
% %%%% integral test of S_minus                          
%          p= -9.747840837847118*1i;                               
%          [q1, q2, n1, n2, alp1, alp2, bet1, bet2] = term_of_gyrotropic_waveguide(EE1, -GG1, HH1, -p);
%          psi = psi1_q__gyrotropic(k_0, k, a_0, EE1, -GG1, HH1, -m, (-p));
%         [B_1, B_2, Cm2, Dm2]  = coefficientsOfContinuousSpectrum(typeOfCylinder, k_0, k, a_0, EE1, -GG1, HH1, MU1, EE, MU, -m, -p, q, psi,  1);
% 
%        p=-p;
%        m=-m;
%        GG1 = -GG1;
%           JM1    = @(r) besselj(m+1, (q1.* r).* k_0);
%           JM2    = @(r) besselj(m+1, (q2.* r).* k_0);
%           Jm1    = @(r) besselj(m, (q1.* r).* k_0);
%           Jm2    = @(r) besselj(m, (q2.* r).* k_0);
%           Jm1_Q1 = @(r) Jm1(r)./ ((q1.* r).* k_0);
%           Jm2_Q2 = @(r) Jm2(r)./ ((q2.* r).* k_0);
%             
%           Ephi_q_minus  = @(r)  B_1.* 1i.* (JM1(r) + (alp1.* m).* Jm1_Q1(r))+...
%                                 B_2.* 1i.* (JM2(r) + (alp2.* m).* Jm2_Q2(r));
%           Hrho_q_minus  = @(r) - 1i * B_1.* (p.*JM1(r) - n1.* (bet1 * m).* Jm1_Q1(r))+...
%                                - 1i * B_2.* (p.*JM2(r) - n2.* (bet2 * m).* Jm2_Q2(r));
%           Hphi_q_minus  = @(r) - B_1.* n1.* (JM1(r) - (bet1 * m).* Jm1_Q1(r))+...
%                                - B_2.* n2.* (JM2(r) - (bet2 * m).* Jm2_Q2(r));
%           Erho_q_minus  = @(r) - B_1.* ((n1.*p + GG1)./EE1.* JM1(r) + (alp1.* m).* Jm1_Q1(r))+...
%                                - B_2.* ((n2.*p + GG1)./EE1.* JM2(r) + (alp2.* m).* Jm2_Q2(r));
%                            
%           
%           A = 1./ (k_0 * (1 - p.^2));
%           A0 = 1./ (k_0 * (1 - p_0.^2));
%           
%           Ez_q0   = @(r) AE_0.* besselj(m0,  k_0.* r* q_0);
%           Hz_q0   = @(r) AH_0.* besselj(m0,  k_0.* r* q_0);
%           dEz_q0  = @(r) AE_0.* k_0.* (q_0) *((besselj(m0,   k_0.* r* q_0) * (m0))./ (k_0.* r* q_0)  - besselj(m0 + 1,   k_0.* r* q_0));
%           dHz_q0  = @(r) AH_0.* k_0.* (q_0) *((besselj(m0,   k_0.* r* q_0) * (m0))./ (k_0.* r* q_0)  - besselj(m0 + 1,   k_0.* r* q_0));
%           
%           Ephi_q_plus0  = @(r)  A0.* ((-p_0.* (m0./r)).* Ez_q0(r) + 1i * dHz_q0(r));
%           Hphi_q_plus0  = @(r)  A0.* (-1i * dEz_q0(r) - (p_0.* (m0./r)).* Hz_q0(r));
%           Hrho_q_plus0  = @(r)  A0.* (((m0./r)).* Ez_q0(r) - 1i * p_0 * dHz_q0(r));
%           Erho_q_plus0  = @(r)  A0.* (-1i * p_0 * dEz_q0(r) - (m0./r).* Hz_q0(r));          
%           
%           S_minus = integral(@(r) (Erho_q_plus0(r).* Hphi_q_minus(r) - Ephi_q_plus0(r).* Hrho_q_minus(r)-...
%               Erho_q_minus(r).* Hphi_q_plus0(r) + Ephi_q_minus(r).* Hrho_q_plus0(r)).* r, 0, a_0);
                                                                               
                                      

        
       aplus2 = (S_plus + S_minus)./ N_p; %%% ????? ?? ????? 1i?
end


    
    
    
    
    
    
    
    
    
