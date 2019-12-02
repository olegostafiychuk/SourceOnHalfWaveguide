function psi2 = psi2_q__gyrotropic(k_0, k, a_0, EE, GG, HH, m, p, q)

%%% w0   - частота поля
%%% a_0  - радиус цилиндра
%%% EE, GG, HH, c, ee - элементы тензора диэлектрической проницаемости,
%%% скорость света в свободном пространстве, диэлектрическая проницаемость
%%% в окружающей среды
%%% AE_0 - коэффициент падения волны с компонентой Ez 
%%% AH_0 - коэффициент падения волны с компонентой Hz 
%%% m    - азимутальный индекс
%%% p   - нормированная постоянная распространения вдоль оси z


%%%%%% вычисляем элементы матрицы рассеяния для внутренней области цилиндра
    mainq = EE.^2 - GG.^2 + EE.* HH - (HH + EE).* p.^2;
    radq = sqrt((HH - EE).^2 * p.^4 + 2 * ((GG.^2).* (HH + EE) - EE.* (HH - EE).^2) * p.^2 +...
        (EE.^2 - GG.^2 - EE.* HH).^2);
    q1 = sqrt(0.5 * (mainq - radq)./ EE);
    q2 = sqrt(0.5 * (mainq + radq)./ EE);

    n1 = -(EE.* (p.* GG).^(-1)).* (p.^2 + q1.^2 + (GG.^2)./ EE - EE);
    n2 = -(EE.* (p * GG).^(-1)).* (p.^2 + q2.^2 + (GG.^2)./ EE - EE);
    
    alp1 = -1 + (p.^2 + q1.^2 - EE)./ GG;
    alp2 = -1 + (p.^2 + q2.^2 - EE)./ GG;
    
%     bet1 = 1 + p./ n1;
%     bet2 = 1 + p./ n2;
    n1bet1 = n1 + p;
    n2bet2 = n2 + p;

    Q1 = (q1.* a_0).* k_0;
    Q2 = (q2.* a_0).* k_0;  
    
    JM1    = besselj(m+1, Q1);
    JM2    = besselj(m+1, Q2);
    Jm1    = besselj(m, Q1);
    Jm2    = besselj(m, Q2);
%     %%%% vary larger JM1 and Jm1
    JM1(abs(imag(Q1))>300) = besselj(m+1, 1i*300);
    Jm1(abs(imag(Q1))>300) = besselj(m,   1i*300);
    Jm1_Q1 = Jm1./ Q1;
    Jm2_Q2 = Jm2./ Q2;
    
    Ez1    = (((1i./HH).* n1).* q1).* Jm1;
    Ez2    = (((1i./HH).* n2).* q2).* Jm2;
    Ephi1  = 1i * (JM1 + (alp1.* m).* Jm1_Q1);
    Ephi2  = 1i * (JM2 + (alp2.* m).* Jm2_Q2);
    Hz1    = - q1.* Jm1;
    Hz2    = - q2.* Jm2;
%     Hphi1  = - n1.* (JM1 - (bet1 * m).* Jm1_Q1);
%     Hphi2  = - n2.* (JM2 - (bet2 * m).* Jm2_Q2);
    Hphi1  = - (n1.* JM1 - (n1bet1 * m).* Jm1_Q1);
    Hphi2  = - (n2.* JM2 - (n2bet2 * m).* Jm2_Q2);
    
    
    %%% компоненты рассеянного поля
%     q = sqrt(1 - p.^2);
%     q = q.* (2*(imag(q) <= 0)-1);
    Q = k.* a_0 * q;
%     Q = Q.* (2*(imag(Q) <= 0)-1); %%% мнимые части S и p должный быть разных знаков
    
    JmmQ  = besselj(m, Q);
    YmmQ  = bessely(m, Q);
    JMmQ  = besselj(m + 1, Q);
    YMmQ  = bessely(m + 1, Q);

%     %%% необходимо, чтобы при стремлении аргумента функции Ханкеля к бесконечности сама функция была ограниченной  
%     H2m  = besselh(m, 2, Q);
%     dH2m = (H2m * m)./ Q  - besselh(m + 1, 2, Q);
    H2m  = JmmQ - 1i * YmmQ;
    dH2m = (H2m * m)./ Q  - (JMmQ - 1i * YMmQ);

    %%% equation for E_rho
    A = 1./ (k_0 * q.^2);   
    Ez_sct  = q.* H2m;
    Hz_sct  = q.* H2m;
    Ephi_sctE  =  q.* A.* (- p.* (m./a_0).* H2m);
    Ephi_sctH  =  q.* A.* (  1i * k * q.* dH2m);
    Hphi_sctE  =  q.* A.* (- 1i * k * q.* dH2m);
    Hphi_sctH  =  q.* A.* (- (p.* (m./a_0)).* H2m);


%%% определитель Delta_1
Delta_1 = (Hz_sct.* (-(Ephi1.* Ez2.* Hphi_sctE) + Ephi1.* Ez_sct.* Hphi2 +...
             Ephi2.* Ez1.* Hphi_sctE - Ephi2.* Ez_sct.* Hphi1 -...
       Ephi_sctE.* Ez1.* Hphi2 +...
             Ephi_sctE.* Ez2.* Hphi1) + Hphi_sctH.* (-(Ephi1.* Ez_sct.* Hz2) +... 
             Ephi2.* Ez_sct.* Hz1 + Ephi_sctE.* Ez1.* Hz2 - Ephi_sctE.* Ez2.* Hz1) +...
        Ephi_sctH.* (-(Ez1.* Hphi_sctE.* Hz2) + Ez2.* Hphi_sctE.* Hz1 +...
             Ez_sct.* Hphi1.* Hz2 - Ez_sct.* Hphi2.* Hz1));
 
%%% определитель Delta_2
%     H1m  = besselh(m, 1, Q);
%     dH1m = (H1m * m)./ Q  - besselh(m + 1, 1, Q);
    H1m  = JmmQ + 1i * YmmQ;
    dH1m = (H1m * m)./ Q  - (JMmQ + 1i * YMmQ);

    %%% equation for E_rho
    A = 1./ (k_0 * q.^2);   
    Ez_sct  = q.* H1m;
    Hz_sct  = q.* H1m;
    Ephi_sctE  =  q.* A.* (- (p.* (m./a_0)).* H1m);
    Ephi_sctH  =  q.* A.* (  1i * k * q.* dH1m);
    Hphi_sctE  =  q.* A.* (- 1i * k * q.* dH1m);
    Hphi_sctH  =  q.* A.* (- (p.* (m./a_0)).* H1m);
    
Delta_2 = (Hz_sct.* (-(Ephi1.* Ez2.* Hphi_sctE) + Ephi1.* Ez_sct.* Hphi2 +...
             Ephi2.* Ez1.* Hphi_sctE - Ephi2.* Ez_sct.* Hphi1 -...
       Ephi_sctE.* Ez1.* Hphi2 +...
             Ephi_sctE.* Ez2.* Hphi1) + Hphi_sctH.* (-(Ephi1.* Ez_sct.* Hz2) +... 
             Ephi2.* Ez_sct.* Hz1 + Ephi_sctE.* Ez1.* Hz2 - Ephi_sctE.* Ez2.* Hz1) +...
        Ephi_sctH.* (-(Ez1.* Hphi_sctE.* Hz2) + Ez2.* Hphi_sctE.* Hz1 +...
             Ez_sct.* Hphi1.* Hz2 - Ez_sct.* Hphi2.* Hz1));
         
%%% определитель Delta_4
    JmQ  = 2 * JmmQ;
    NmQ  = 2 * 1i * YmmQ;

    %%% equation for E_rho
    A = 1./ (k_0 * q.^2);   
    Ez_sct  = q.* JmQ;
    Hz_sct  = q.* NmQ;
    Ephi_sctE  =  q.* A.* (- (p.* (m./a_0)).* JmQ);
    Hphi_sctE  =  q.* A.* (- 1i * k * q.* ((JmQ * m)./ Q  - 2 * JMmQ));
    Ephi_sctH  =  q.* A.* (  1i * k * q.* ((NmQ * m)./ Q  - 2 * 1i * YMmQ));
    Hphi_sctH  =  q.* A.* (-  p.* (m./a_0).* NmQ);
    
Delta_4 = (Hz_sct.* (-(Ephi1.* Ez2.* Hphi_sctE) + Ephi1.* Ez_sct.* Hphi2 +...
             Ephi2.* Ez1.* Hphi_sctE - Ephi2.* Ez_sct.* Hphi1 -...
       Ephi_sctE.* Ez1.* Hphi2 +...
             Ephi_sctE.* Ez2.* Hphi1) + Hphi_sctH.* (-(Ephi1.* Ez_sct.* Hz2) +... 
             Ephi2.* Ez_sct.* Hz1 + Ephi_sctE.* Ez1.* Hz2 - Ephi_sctE.* Ez2.* Hz1) +...
        Ephi_sctH.* (-(Ez1.* Hphi_sctE.* Hz2) + Ez2.* Hphi_sctE.* Hz1 +...
             Ez_sct.* Hphi1.* Hz2 - Ez_sct.* Hphi2.* Hz1));

sqrtDelt = sqrt((Delta_4 + Delta_1 - Delta_2).^2 + 4 * Delta_2.* Delta_1);
% sqrtDelt = (real(sqrtDelt) + 1i* abs(imag(sqrtDelt))).* (2*(q <= 1)-1).*...
%            (sign(real(p))-sign(imag(p))).* (2*(GG <= 0)-1);
sqrtDelt = (real(sqrtDelt) + 1i* abs(imag(sqrtDelt)));
sqrtDelt = sqrtDelt.* (2*(q <= 1)-1).*...
    (sign(real(p))-sign(imag(p))).* (2*(GG <= 0)-1);
psi2 = (-(Delta_4 + Delta_1 - Delta_2)-...
    sqrtDelt)./ (2 * Delta_2);

         
% B1 = (Hz_sct.* (Ephi2.* Ez_inc.* Hphi_sctE - Ephi2.* Ez_sct.* Hphi_inc -...
%              Ephi_inc.* Ez2.* Hphi_sctE + Ephi_inc.* Ez_sct.* Hphi2 +...
%              Ephi_sctE.* Ez2.* Hphi_inc - Ephi_sctE.* Ez_inc.* Hphi2) +...
%         Hphi_sctH.* (Ephi2.* Ez_sct.* Hz_inc - Ephi_inc.* Ez_sct.* Hz2 -...
%              Ephi_sctE.* Ez2.* Hz_inc + Ephi_sctE.* Ez_inc.* Hz2) +...
%         Ephi_sctH.* (Ez2.* Hphi_sctE.* Hz_inc - Ez_inc.* Hphi_sctE.* Hz2 -...
%              Ez_sct.* Hphi2.* Hz_inc + Ez_sct.* Hphi_inc.* Hz2))./Delta;
%     
%          
% B2 = (Hz_sct.* (-(Ephi1.* Ez_inc.* Hphi_sctE) + Ephi1.* Ez_sct.* Hphi_inc +... 
%              Ephi_inc.* Ez1.* Hphi_sctE - Ephi_inc.* Ez_sct.* Hphi1 -...
%              Ephi_sctE.* Ez1.* Hphi_inc + Ephi_sctE.* Ez_inc.* Hphi1) +...
%         Hphi_sctH.* (-(Ephi1.* Ez_sct.* Hz_inc) + Ephi_inc.* Ez_sct.* Hz1 +...
%              Ephi_sctE.* Ez1.* Hz_inc - Ephi_sctE.* Ez_inc.* Hz1) +...
%         Ephi_sctH.* (-(Ez1.* Hphi_sctE.* Hz_inc) + Ez_inc.* Hphi_sctE.* Hz1 +...
%              Ez_sct.* Hphi1.* Hz_inc - Ez_sct.* Hphi_inc.* Hz1))./ Delta; 
%              
%          
% DE = (Hz_sct.* (Ephi1.* Ez2.* Hphi_inc - Ephi1.* Ez_inc.* Hphi2 -...
%              Ephi2.* Ez1.* Hphi_inc + Ephi2.* Ez_inc.* Hphi1 +...
%        Ephi_inc.* Ez1.* Hphi2 -...
%              Ephi_inc.* Ez2.* Hphi1) + Hphi_sctH.* (-(Ephi1.* Ez2.* Hz_inc) +...
%              Ephi1.* Ez_inc.* Hz2 + Ephi2.* Ez1.* Hz_inc - Ephi2.* Ez_inc.* Hz1 -...
%              Ephi_inc.* Ez1.* Hz2 + Ephi_inc.* Ez2.* Hz1) +...
%         Ephi_sctH.* (-(Ez1.* Hphi2.* Hz_inc) + Ez1.* Hphi_inc.* Hz2 +...
%              Ez2.* Hphi1.* Hz_inc - Ez2.* Hphi_inc.* Hz1 - Ez_inc.* Hphi1.* Hz2 +...
%              Ez_inc.* Hphi2.* Hz1))./Delta; 
%              
%          
% DH = (Hz_inc.* (Ephi1.* Ez2.* Hphi_sctE - Ephi1.* Ez_sct.* Hphi2 -...
%              Ephi2.* Ez1.* Hphi_sctE + Ephi2.* Ez_sct.* Hphi1) -...
%         (Ephi2.* Hz1 - Ephi1.* Hz2).* (Ez_sct.* Hphi_inc - Ez_inc.* Hphi_sctE) +...
%         Ephi_inc.* (Ez1.* Hphi_sctE.* Hz2 - Ez2.* Hphi_sctE.* Hz1 -...
%        Ez_sct.* Hphi1.* Hz2 +...
%              Ez_sct.* Hphi2.* Hz1) + Ephi_sctE.* (Ez1.* Hphi2.* Hz_inc -...
%              Ez1.* Hphi_inc.* Hz2 - Ez2.* Hphi1.* Hz_inc + Ez2.* Hphi_inc.* Hz1 +...
%              Ez_inc.* Hphi1.* Hz2 - Ez_inc.* Hphi2.* Hz1))./Delta;

         
%          DH = Q1;
