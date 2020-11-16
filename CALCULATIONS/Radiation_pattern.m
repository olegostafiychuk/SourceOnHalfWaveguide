clearvars
clc

% вычисление диаграммы направленности несимметричного источника
typeOfCylinder = 'Gyrotropic';
% typeOfCylinder = 'Isotropic';
systemParameters

% EE1 =  1.00001;
% GG1 = 0.000001;
% HH1 =  1.00001;

step = 1/5000;
teta1 = (0.0001:step:0.49999)*pi;
teta2 = (0.50001:step:0.9999)*pi;
%teta2 = teta1;
q1 = sin(teta1);
q2 = sin(teta2);
p1 = sqrt(1 - q1.^2);%cos(teta1);%
p2 = sqrt(1 - q2.^2);%cos(teta2);%
z = 0;
psi1   = psi1_q__gyrotropic(k_0, k_0, a_0, EE1,  GG1, HH1,  m, p1, q1);
psi2   = psi2_q__gyrotropic(k_0, k_0, a_0, EE1,  GG1, HH1,  m, p1, q1);
P_plus_sum = zeros(1,size(teta1,2));
P_minus_sum = zeros(1,size(teta2,2));
% P_plus_sum = zeros(1,size(q1,2));
% P_minus_sum = zeros(1,size(q2,2));
for i = 1:2
    for k =1:2
        switch k
            case 1
                a_plus_bet = a_sma_of_continuousWaves_alpha1(typeOfCylinder, q1, p_0, p1, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, z, j_f, j_z, d);
                a_minus_bet = a_sma_of_continuousWaves_alpha1(typeOfCylinder, q2, p_0, -p2, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, z, j_f, j_z, d);
                [B_1, B_2, Cm_2_plus_bet, Dm_2_plus_bet]  = coefficientsOfContinuousSpectrum(typeOfCylinder, k_0, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, p1, q1, psi1,  1);
                psi_bet = psi1_q__gyrotropic(k_0, k_0, a_0, EE1,  GG1, HH1,  m, -p2, q2);
                [B_1, B_2, Cm_2_minus_bet, Dm_minus_2_bet]  = coefficientsOfContinuousSpectrum(typeOfCylinder, k_0, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, -p2, q2, psi_bet,  1);
                
            case 2
                a_plus_bet = a_sma_of_continuousWaves_alpha2(typeOfCylinder, q1, p_0, p1, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, z, j_f, j_z, d);
                a_minus_bet = a_sma_of_continuousWaves_alpha2(typeOfCylinder, q2, p_0, -p2, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, z, j_f, j_z, d);
                [B_1, B_2, Cm_2_plus_bet, Dm_2_plus_bet]  = coefficientsOfContinuousSpectrum_forPsi2(typeOfCylinder, k_0, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, p1, q1, psi2,  1);
                psi_bet = psi2_q__gyrotropic(k_0, k_0, a_0, EE1,  GG1, HH1,  m, -p2, q2);
                [B_1, B_2, Cm_2_minus_bet, Dm_minus_2_bet]  = coefficientsOfContinuousSpectrum_forPsi2(typeOfCylinder, k_0, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, -p2, q2, psi_bet,  1);
                
        end
        switch i
            case 1
                a_plus_alph = a_sma_of_continuousWaves_alpha1(typeOfCylinder, q1, p_0, p1, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, z, j_f, j_z, d);
                a_minus_alph = a_sma_of_continuousWaves_alpha1(typeOfCylinder, q2, p_0, -p2, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, z, j_f, j_z, d);
                [B_1, B_2, Cm_2_plus_alph, Dm_2_plus_alph]  = coefficientsOfContinuousSpectrum(typeOfCylinder, k_0, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, p1, q1, psi1,  1);
                psi_alph = psi1_q__gyrotropic(k_0, k_0, a_0, EE1,  GG1, HH1,  m, -p2, q2);
                [B_1, B_2, Cm_2_minus_alph, Dm_2_minus_alph]  = coefficientsOfContinuousSpectrum(typeOfCylinder, k_0, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, -p2, q2, psi_alph,  1);
            case 2
                a_plus_alph = a_sma_of_continuousWaves_alpha2(typeOfCylinder, q1, p_0, p1, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, z, j_f, j_z, d);
                a_minus_alph = a_sma_of_continuousWaves_alpha2(typeOfCylinder, q2, p_0, -p2, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, z, j_f, j_z, d);
                [B_1, B_2, Cm_2_plus_alph, Dm_2_plus_alph]  = coefficientsOfContinuousSpectrum_forPsi2(typeOfCylinder, k_0, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, p1, q1, psi2,  1);
                psi_alph = psi2_q__gyrotropic(k_0, k_0, a_0, EE1,  GG1, HH1,  m, -p2, q2);
                [B_1, B_2, Cm_2_minus_alph, Dm_2_minus_alph]  = coefficientsOfContinuousSpectrum_forPsi2(typeOfCylinder, k_0, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, -p2, q2, psi_alph,  1);
                
        end
        
        %if(k==i && i==2)
        P_plus_sum = P_plus_sum + a_plus_alph.*conj(a_plus_bet).*...
           (Cm_2_plus_alph.*conj(Cm_2_plus_bet) + Dm_2_plus_alph.*conj(Dm_2_plus_bet));
%         P_plus_sum = P_plus_sum + Cm_2_plus_alph.*conj(Cm_2_plus_bet) + Dm_2_plus_alph.*conj(Dm_2_plus_bet);
       
       
        P_minus_sum = P_minus_sum + a_minus_alph.*conj(a_minus_bet).*...
            (Cm_2_minus_alph.*conj(Cm_2_minus_bet) + Dm_2_minus_alph.*conj(Dm_minus_2_bet));%.*psi_alph.*conj(psi_bet);
%         P_minus_sum = P_minus_sum + Cm_2_minus_alph.*conj(Cm_2_minus_bet) + Dm_2_minus_alph.*conj(Dm_minus_2_bet);
       % end
        
    end
end

D_plus = c/(2*pi*k_0^2).*(cot(teta1)).^2.*P_plus_sum; % (пока только для
%положительного направления оси z) %
D_minus = c/(2*pi*k_0^2).*(cot(teta2)).^2.*P_minus_sum; % для отрицательного направления оси z
D_max = max(max(abs(D_plus)),max(abs(D_minus)));

% figure(1)
% %subplot(1,2,1)
% plot(teta1/pi,real(D_plus/D_max), teta1/pi, imag(D_plus/D_max),'--');  hold on;
% xlabel('teta (\pi rad)'); ylabel('S_r/max(|S_r|)');
% %subplot(1,2,2)
% plot(teta2/pi,real(D_minus/D_max), teta2/pi, imag(D_minus/D_max),'--'); hold on;
% %xlabel('teta')

teta = [teta1  teta2];
D1 = [D_plus zeros(1,size(teta2,2))];
D2 = [zeros(1,size(teta1,2)) D_minus];
D = (D1 + D2)/D_max;
figure(1)
plot(teta/pi, real(D), 'k'); hold on;

% figure(8)
% plot(teta1, D_plus./D_max); hold on
% teta_inv = pi.*ones(1,length(teta2))-teta2;
% plot(teta_inv, D_minus./D_max);

xlabel('teta (\pi rad)'); ylabel('S_r/max(|S_r|)');

% Построение в полярных координатах
teta = [teta1  teta2];
D1 = [D_plus zeros(1,size(teta2,2))];
D2 = [zeros(1,size(teta1,2)) D_minus];
D = (D1 + D2)/D_max;
figure(2)
polar(teta, D,'k'); hold on
polar(-teta, D,'k');

% polar(teta1,real(D_plus/D_max));  %hold on;
% polar(teta2,real(D_minus/D_max));

















