function [B1,B2,Cm2, Dm] = coefficientsOfDiscreteSpectrum(typeOfCylinder, k_0, k, a_0, EE1, GG1, HH1, MU1, EE, MU, m, p, q, psi, Dm2)
%%% вычислеяются коэффициенты собственных волн цилиндрического волновода
%%% (как дискретного, так и непрерывного спектра) от коэффициента Dm2
%%% (принят в качестве известного)

switch(typeOfCylinder)
    case 'PerfectConductivity'
        q1 = sqrt(EE1*MU1/(EE*MU) -  p.^2);
        Q1 = (q1.* a_0).* k_0;
        q = sqrt(1 - p.^2);
        %     q = real(q) - 1i * abs(imag(q));
        q = q.* (2*(imag(q) <= 0)-1);
        Q = k_0.* a_0 * q;
        %     Q = Q.* (2*(imag(Q) <= 0)-1); %%% мнимые части S и p должный быть разных знаков
        %     %%% необходимо, чтобы при стремлении аргумента функции Ханкеля к бесконечности сама функция была ограниченной  
            H1m  = 0;
            dH1m = 0;
            
            H2m  = besselh(m, 2, Q);
            dH2m = (H2m * m)./ Q  - besselh(m + 1, 2, Q);
            
            Jm1    = besselj(m, Q1);
            dJm1    = (m./ Q1).* Jm1 - besselj(m+1, Q1); 
            
            A = 1./ (k_0 * EE*MU * (1 - p.^2)); 
            
            A_1 = 1./ (k_0 * (EE1*MU1 - EE*MU * p.^2));
            M11 = Jm1;
            M12 = 0;
            
            M13 = -(psi.* H1m + H2m);
            M21 = 0;
            M22 = Jm1;
            M23 = 0;
            %%%% Ephi
            M31 = -A_1.* ((p * m./ a_0).* Jm1);
            M32 =  A_1.* ( 1i * MU1 * k_0 * q1.* dJm1);
            M33 =  A.* ((p.* (m./a_0)).* (psi.* H1m + H2m));
            %     %%%% Hphi
            %     M31 =  A_1.* (-1i * EE1 * k_0 * q1.* dJm1);
            %     M32 =  A_1.* ((-p * m./ a_0).* Jm1);
            %     M33 =  A.* ( 1i * EE * k_0 * q.* (psi.* dH1m + dH2m));
            
            
            %%% правая часть
            R1  = 0;
            R2  = -psi.* H1m + H2m;
            R3  =  A.* (1i * MU * k_0 * q.* (-psi.* dH1m + dH2m)); %%%% Ephi
            %     R3  =  A.* ((-p.* (m./a_0)).*(-psi.* H1m + H2m)); %%%% Hphi
            
            %%% Вычислеям коэффициенты 
            %%% Dm равно детерминанту системы
            Dm = -M13.* M22.* M31 + M12.* M23.* M31 + M13.* M21.* M32 -...
                M11.* M23.* M32 - M12.* M21.* M33 + M11.* M22.* M33;
            
            B1 = (-M13.* M22.* R3 + M12.* M23.* R3 + M13.* R2.* M32 -...
                R1.* M23.* M32 - M12.* R2.* M33 + R1.* M22.* M33);
            
            B2 = (-M13.* R2.* M31 + R1.* M23.* M31 + M13.* M21.* R3 -...
                M11.* M23.* R3 - R1.* M21.* M33 + M11.* R2.* M33);
            
            Cm2 = (-R1.* M22.* M31 + M12.* R2.* M31 + R1.* M21.* M32 -...
                M11.* R2.* M32 - M12.* M21.* R3 + M11.* M22.* R3);
            
            
            B1 = B1./Dm;
            B2 = B2./Dm;
            Cm2= Cm2./Dm;
            Dm = 1;
    case 'Isotropic'
        q1 = sqrt(EE1*MU1/(EE*MU) -  p.^2);
        Q1 = (q1.* a_0).* k_0;
%         q = sqrt(1 - p.^2);
        %     q = real(q) - 1i * abs(imag(q));
            q = q.* (2*(imag(q) <= 0)-1);
            Q = k_0.* a_0 * q;
            %     Q = Q.* (2*(imag(Q) <= 0)-1); %%% мнимые части S и p должный быть разных знаков
            
            H1m  = 0;
            dH1m = 0;
            
            H2m  = besselh(m, 2, Q);
            dH2m = (H2m * m)./ Q  - besselh(m + 1, 2, Q);
            
            Jm1    = besselj(m, Q1);
            dJm1    = (m./ Q1).* Jm1 - besselj(m+1, Q1); 
            
            A = 1./ (k_0 * EE*MU * q.^2); 
            
            A_1 = 1./ (k_0 * (EE1*MU1 - EE*MU * p.^2));
            M11 = q1.* Jm1;
            M12 = 0;            
            M13 = -q.* (psi.* H1m + H2m);
            M21 = 0;
            M22 = q1.* Jm1;
            M23 = 0;
            %%%% Ephi
            M31 = -A_1.* q1.* ((p * m./ a_0).* Jm1);
            M32 =  A_1.* q1.* ( 1i * MU1 * k_0 * q1.* dJm1);
            M33 =  A.* q.* ((p.* (m./a_0)).* (psi.* H1m + H2m));
            %     %%%% Hphi
            %     M31 =  A_1.* (-1i * EE1 * k_0 * q1.* dJm1);
            %     M32 =  A_1.* ((-p * m./ a_0).* Jm1);
            %     M33 =  A.* ( 1i * EE * k_0 * q.* (psi.* dH1m + dH2m));
            
            
            %%% правая часть
            R1  = 0;
            R2  = q.* (-psi.* H1m + H2m);
            R3  =  q.* A.* (1i * MU * k_0 * q.* (-psi.* dH1m + dH2m)); %%%% Ephi
            %     R3  =  A.* ((-p.* (m./a_0)).*(-psi.* H1m + H2m)); %%%% Hphi
            
            %%% Вычислеям коэффициенты 
            %%% Dm равно детерминанту системы
            Dm = -M13.* M22.* M31 + M12.* M23.* M31 + M13.* M21.* M32 -...
                M11.* M23.* M32 - M12.* M21.* M33 + M11.* M22.* M33;
            
            B1 = (-M13.* M22.* R3 + M12.* M23.* R3 + M13.* R2.* M32 -...
                R1.* M23.* M32 - M12.* R2.* M33 + R1.* M22.* M33);
            
            B2 = (-M13.* R2.* M31 + R1.* M23.* M31 + M13.* M21.* R3 -...
                M11.* M23.* R3 - R1.* M21.* M33 + M11.* R2.* M33);
            
            Cm2 = (-R1.* M22.* M31 + M12.* R2.* M31 + R1.* M21.* M32 -...
                M11.* R2.* M32 - M12.* M21.* R3 + M11.* M22.* R3);
            
            
            B1 = B1./Dm;
            B2 = B2./Dm;
            Cm2= Cm2./Dm;
            Dm = Dm./ Dm;
            
    case 'Gyrotropic'
%         q = sqrt(1 - p.^2);
        %     q = real(q) - 1i * abs(imag(q));
            q = q.* (2*(imag(q) <= 0)-1);
            Q = k_0.* a_0 * q;
            %     Q = Q.* (2*(imag(Q) <= 0)-1); %%% мнимые части S и p должный быть разных знаков
            %     %%% необходимо, чтобы при стремлении аргумента функции Ханкеля к бесконечности сама функция была ограниченной  

            H1m  = 0;
            dH1m = 0;
            
            H2m  = besselh(m, 2, Q);
            dH2m = (H2m * m)./ Q  - besselh(m + 1, 2, Q);
            
            mainq = EE1.^2 - GG1.^2 + EE1.* HH1 - (HH1 + EE1).* p.^2;
            radq = sqrt((HH1 - EE1).^2 * p.^4 + 2 * ((GG1.^2).* (HH1 + EE1) - EE1.* (HH1 - EE1).^2) * p.^2 +...
                (EE1.^2 - GG1.^2 - EE1.* HH1).^2);
            q1 = sqrt(0.5 * (mainq - radq)./ EE1);
            q2 = sqrt(0.5 * (mainq + radq)./ EE1);
            
            n1 = -(EE1.* (p.* GG1).^(-1)).* (p.^2 + q1.^2 + (GG1.^2)./ EE1 - EE1);
            n2 = -(EE1.* (p * GG1).^(-1)).* (p.^2 + q2.^2 + (GG1.^2)./ EE1 - EE1);
            
            alp1 = -1 + (p.^2 + q1.^2 - EE1)./ GG1;
            alp2 = -1 + (p.^2 + q2.^2 - EE1)./ GG1;
            
%           bet1 = 1 + p./ n1;
%           bet2 = 1 + p./ n2;
            n1bet1 = n1 + p;
            n2bet2 = n2 + p;
            
            Q1 = (q1.* a_0).* k_0;
            Q2 = (q2.* a_0).* k_0;
            
            JM1    = besselj(m+1, Q1);
            JM2    = besselj(m+1, Q2);
            Jm1    = besselj(m, Q1);
            Jm2    = besselj(m, Q2);
%             %%%% vary larger JM1 and Jm1
            JM1(abs(imag(Q1))>300) = besselj(m+1, 1i*300);
            Jm1(abs(imag(Q1))>300) = besselj(m, 1i*300);
            
            Jm1_Q1 = Jm1./ Q1;
            Jm2_Q2 = Jm2./ Q2;
            
            Ez1    = (((1i./HH1).* n1).* q1).* Jm1;
            Ez2    = (((1i./HH1).* n2).* q2).* Jm2;
            Ephi1  = 1i * (JM1 + (alp1.* m).* Jm1_Q1);
            Ephi2  = 1i * (JM2 + (alp2.* m).* Jm2_Q2);
            Hz1    = - q1.* Jm1;
            Hz2    = - q2.* Jm2;
% %             Hphi1  = - n1.* (JM1 - (bet1 * m).* Jm1_Q1);
% %             Hphi2  = - n2.* (JM2 - (bet2 * m).* Jm2_Q2);
%             Hphi1  = - (n1.* JM1 - (n1bet1 * m).* Jm1_Q1);
%             Hphi2  = - (n2.* JM2 - (n2bet2 * m).* Jm2_Q2)
            A = 1./ (k_0 * q.^2); 
            
            M11 = Ez1;
            M12 = Ez2;
            M13 = -q.* (psi.* H1m + H2m);
            
            M21 = Hz1;
            M22 = Hz2;
            M23 = 0;
            %%%% Ephi
            M31 = Ephi1;
            M32 = Ephi2;
            M33 =  q.* A.* ((p.* (m./a_0)).* (psi.* H1m + H2m));
            %     %%%% Hphi
            %     M31 =  A_1.* (-1i * EE1 * k_0 * q1.* dJm1);
            %     M32 =  A_1.* ((-p * m./ a_0).* Jm1);
            %     M33 =  A.* ( 1i * EE * k_0 * q.* (psi.* dH1m + dH2m));
            
            
            %%% правая часть
            R1  = 0;
            R2  = q.* (-psi.* H1m + H2m);
            R3  = q.* A.* (1i * MU * k_0 * q.* (-psi.* dH1m + dH2m)); %%%% Ephi
            %     R3  =  A.* ((-p.* (m./a_0)).*(-psi.* H1m + H2m)); %%%% Hphi
            
            %%% Вычислеям коэффициенты 
            %%% Dm равно детерминанту системы
            Dm = -M13.* M22.* M31 + M12.* M23.* M31 + M13.* M21.* M32 -...
                M11.* M23.* M32 - M12.* M21.* M33 + M11.* M22.* M33;
            
            B1 = (-M13.* M22.* R3 + M12.* M23.* R3 + M13.* R2.* M32 -...
                R1.* M23.* M32 - M12.* R2.* M33 + R1.* M22.* M33);
            %%%% for vary larger JM1 and Jm1
            B1(abs(imag(Q1))>300) = 0;
            
            B2 = (-M13.* R2.* M31 + R1.* M23.* M31 + M13.* M21.* R3 -...
                M11.* M23.* R3 - R1.* M21.* M33 + M11.* R2.* M33);
            
            Cm2 = (-R1.* M22.* M31 + M12.* R2.* M31 + R1.* M21.* M32 -...
                M11.* R2.* M32 - M12.* M21.* R3 + M11.* M22.* R3);
            
            
            if(Dm ~= 0)
                B1 = B1./ Dm;
                B2 = B2./ Dm;
                Cm2= Cm2./Dm;
                Dm = Dm./ Dm;
            end

% %             B1 = B1./ Cm2;
% %             B2 = B2./ Cm2;
% %             Cm2= Cm2./Cm2;
% %             Dm = Dm./ Cm2;
end
    
% %%%%%%%%%% решаем неоднородную систему линейных уравнений
% %%%% левая (однородная) часть системы
% mat_boundConditions = [ M11      M12      M13 ;
%                         M21      M22      M23 ;
%                         M31      M32      M33 ];
% %%%% правая (неоднородная) часть системы
% rightSide = [R1;
%              R2; 
%              R3].* Dm2;
% A = mat_boundConditions \ rightSide;
% B1  = A(1);
% B2  = A(2);
% Cm2 = A(3);


