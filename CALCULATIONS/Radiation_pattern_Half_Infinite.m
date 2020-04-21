clearvars
clc

% вычисление диаграммы направленности несимметричного источника
typeOfCylinder = 'Gyrotropic';
% typeOfCylinder = 'Isotropic';
systemParameters

% EE1 =  1.00001;
% GG1 = 0.000001;
% HH1 =  1.00001;
% waveguideParameters.EEinner = 1.00001;
% waveguideParameters.GGinner = 0.000001;
% waveguideParameters.HHinner = 1.00001;

% saveVarsMat = load('C:\Users\Олег\Dropbox\PAPERS\URSI GASS 2020\coefficients_0995wUH_n3e12_L_10d.mat');
saveVarsMat = load('C:\Users\Олег\Dropbox\PAPERS\URSI GASS 2020\coefficients_0025wH_n1e13_L_2lambda0.mat');
p = saveVarsMat.p_cs; 
q = saveVarsMat.q_cs;
a_cs1_minus = saveVarsMat.a_cs_alfa1_back_prop;
a_cs2_minus = saveVarsMat.a_cs_alfa2_back_prop;
b_cs1_minus = saveVarsMat.b_p_field_1_back_prop;
b_cs2_minus = saveVarsMat.b_p_field_2_back_prop;
c_cs1 = saveVarsMat.a_p_Ebeam_forw_prop;
c_cs2 = saveVarsMat.a_p_Hbeam_forw_prop;
L = saveVarsMat.LL;
clear saveVarsMat;
ind = find(q < 1);
ind = ind(end);
q = q(2:ind);
p = p(2:ind);
a_cs1_minus = a_cs1_minus(2:ind);
a_cs2_minus = a_cs2_minus(2:ind);
b_cs1_minus = b_cs1_minus(2:ind);
b_cs2_minus = b_cs2_minus(2:ind);
c_cs1 = c_cs1(2:ind);
c_cs2 = c_cs2(2:ind);

coeffsOfField = coeffsOfField_of_Continues_discreteRepr(typeOfCylinder, q, p, waveguideParameters, sourceParameters);

psi_forward_1 = coeffsOfField.psi_forward_1;
psi_backward_1 = coeffsOfField.psi_backward_1;
psi_forward_1_transp = coeffsOfField.psi_forward_1_transp;
psi_backward_1_transp = coeffsOfField.psi_backward_1_transp;
psi_forward_2 = coeffsOfField.psi_forward_2;
psi_backward_2 = coeffsOfField.psi_backward_2;
psi_forward_2_transp = coeffsOfField.psi_forward_2_transp;
psi_backward_2_transp = coeffsOfField.psi_backward_2_transp;
B_1_forward_1 = coeffsOfField.B_1_forward_1;
B_2_forward_1 = coeffsOfField.B_2_forward_1;
Cm2_forward_1 = coeffsOfField.Cm2_forward_1;
Dm2_forward_1 = coeffsOfField.Dm2_forward_1;
B_1_backward_1 = coeffsOfField.B_1_backward_1;
B_2_backward_1 = coeffsOfField.B_2_backward_1;
Cm2_backward_1 = coeffsOfField.Cm2_backward_1;
Dm2_backward_1 = coeffsOfField.Dm2_backward_1;
B_1_forward_2 = coeffsOfField.B_1_forward_2;
B_2_forward_2 = coeffsOfField.B_2_forward_2;
Cm2_forward_2 = coeffsOfField.Cm2_forward_2;
Dm2_forward_2 = coeffsOfField.Dm2_forward_2;
B_1_backward_2 = coeffsOfField.B_1_backward_2;
B_2_backward_2 = coeffsOfField.B_2_backward_2;
Cm2_backward_2 = coeffsOfField.Cm2_backward_2;
Dm2_backward_2 = coeffsOfField.Dm2_backward_2;


%%%%%%%%%%%%%%%%%%%%

teta1 = asin(q);
teta2 = pi/2 + asin(q);
q1 = sin(teta1);
q2 = sin(teta2);
p1 = cos(teta1);%sqrt(1 - q1.^2);%cos(teta1);
p2 = cos(teta2);%sqrt(1 - q2.^2);%
% plot(q2, q); hold on;
% plot(teta2, p2)

% a_cs1_minus = a_cs1_minus(end:-1:1);
% a_cs2_minus = a_cs2_minus(end:-1:1);
% b_cs1_minus = b_cs1_minus(end:-1:1);
% b_cs2_minus = b_cs2_minus(end:-1:1);
% psi_backward_1 = psi_backward_1(end:-1:1);
% psi_backward_2 = psi_backward_2(end:-1:1);
% Cm2_backward_1 = Cm2_backward_1(end:-1:1);
% Cm2_backward_2 = Cm2_backward_2(end:-1:1);
% Dm2_backward_1 = Dm2_backward_1(end:-1:1);
% Dm2_backward_2 = Dm2_backward_2(end:-1:1);

P_plus_sum = zeros(1,size(q,2));
P_minus_sum = zeros(1,size(q,2));
for i = 1:2
    for k =1:2
        switch k
            case 1
                %a_plus_bet = a_sma_of_continuousWaves_alpha1(typeOfCylinder, q1, p_0, p1, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, z, j_f, j_z, d);
%                 a_minus_bet = a_sma_of_continuousWaves_alpha1(typeOfCylinder, q2, p_0, -p2, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, z, j_f, j_z, d);
                  %[B_1, B_2, Cm_2_plus_bet, Dm_2_plus_bet]  = coefficientsOfContinuousSpectrum(typeOfCylinder, k_0, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, p1, q1, psi1,  1);
%                 psi_bet = psi1_q__gyrotropic(k_0, k_0, a_0, EE1,  GG1, HH1,  m, -p2, q2);
%                 [B_1, B_2, Cm_2_minus_bet, Dm_minus_2_bet]  = coefficientsOfContinuousSpectrum(typeOfCylinder, k_0, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, -p2, q2, psi_bet,  1);
                                
                a_minus_bet = a_cs1_minus;
                b_minus_bet = b_cs1_minus;
                psi_bet = psi_backward_1;
                Cm_2_minus_bet = Cm2_backward_1;
                Dm_minus_2_bet = Dm2_backward_1;
%                 psi_bet = psi1_q__gyrotropic(k_0, k_0, a_0, EE1,  GG1, HH1,  m, p2, q2);
%                 [B_1, B_2, Cm_2_minus_bet, Dm_minus_2_bet]  = coefficientsOfContinuousSpectrum(typeOfCylinder, k_0, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, p2, q2, psi_bet,  1);
                                
              
            case 2
                %a_plus_bet = a_sma_of_continuousWaves_alpha2(typeOfCylinder, q1, p_0, p1, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, z, j_f, j_z, d);
%                 a_minus_bet = a_sma_of_continuousWaves_alpha2(typeOfCylinder, q2, p_0, -p2, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, z, j_f, j_z, d);
                %[B_1, B_2, Cm_2_plus_bet, Dm_2_plus_bet]  = coefficientsOfContinuousSpectrum_forPsi2(typeOfCylinder, k_0, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, p1, q1, psi2,  1);
%                 psi_bet = psi2_q__gyrotropic(k_0, k_0, a_0, EE1,  GG1, HH1,  m, -p2, q2);
%                 [B_1, B_2, Cm_2_minus_bet, Dm_minus_2_bet]  = coefficientsOfContinuousSpectrum_forPsi2(typeOfCylinder, k_0, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, -p2, q2, psi_bet,  1);
                
                a_minus_bet = a_cs2_minus;
                b_minus_bet = b_cs2_minus;
                psi_bet = psi_backward_2;
                Cm_2_minus_bet = Cm2_backward_2;
                Dm_minus_2_bet = Dm2_backward_2;
%                 psi_bet = psi2_q__gyrotropic(k_0, k_0, a_0, EE1,  GG1, HH1,  m, p2, q2);
%                 [B_1, B_2, Cm_2_minus_bet, Dm_minus_2_bet]  = coefficientsOfContinuousSpectrum_forPsi2(typeOfCylinder, k_0, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, p2, q2, psi_bet,  1);
                
                
        end
        switch i
            case 1
                %a_plus_alph = a_sma_of_continuousWaves_alpha1(typeOfCylinder, q1, p_0, p1, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, z, j_f, j_z, d);
%                 a_minus_alph = a_sma_of_continuousWaves_alpha1(typeOfCylinder, q2, p_0, -p2, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, z, j_f, j_z, d);
                %[B_1, B_2, Cm_2_plus_alph, Dm_2_plus_alph]  = coefficientsOfContinuousSpectrum(typeOfCylinder, k_0, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, p1, q1, psi1,  1);
%                 psi_alph = psi1_q__gyrotropic(k_0, k_0, a_0, EE1,  GG1, HH1,  m, -p2, q2);
%                 [B_1, B_2, Cm_2_minus_alph, Dm_2_minus_alph]  = coefficientsOfContinuousSpectrum(typeOfCylinder, k_0, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, -p2, q2, psi_alph,  1);
            
                a_minus_alph = a_cs1_minus;
                b_minus_alph = b_cs1_minus;
                psi_alph = psi_backward_1;
                Cm_2_minus_alph = Cm2_backward_1;
                Dm_minus_2_alph = Dm2_backward_1;
%                 psi_alph = psi1_q__gyrotropic(k_0, k_0, a_0, EE1,  GG1, HH1,  m, p2, q2);
%                 [B_1, B_2, Cm_2_minus_alph, Dm_2_minus_alph]  = coefficientsOfContinuousSpectrum(typeOfCylinder, k_0, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, p2, q2, psi_alph,  1);
            
            
            case 2
                %a_plus_alph = a_sma_of_continuousWaves_alpha2(typeOfCylinder, q1, p_0, p1, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, z, j_f, j_z, d);
%                 a_minus_alph = a_sma_of_continuousWaves_alpha2(typeOfCylinder, q2, p_0, -p2, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, z, j_f, j_z, d);
                %[B_1, B_2, Cm_2_plus_alph, Dm_2_plus_alph]  = coefficientsOfContinuousSpectrum_forPsi2(typeOfCylinder, k_0, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, p1, q1, psi2,  1);
%                 psi_alph = psi2_q__gyrotropic(k_0, k_0, a_0, EE1,  GG1, HH1,  m, -p2, q2);
%                 [B_1, B_2, Cm_2_minus_alph, Dm_2_minus_alph]  = coefficientsOfContinuousSpectrum_forPsi2(typeOfCylinder, k_0, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, -p2, q2, psi_alph,  1);
%                 
                a_minus_alph = a_cs2_minus;
                b_minus_alph = b_cs2_minus;
                psi_alph = psi_backward_2;
                Cm_2_minus_alph = Cm2_backward_2;
                Dm_minus_2_alph = Dm2_backward_2;
%                 psi_alph = psi2_q__gyrotropic(k_0, k_0, a_0, EE1,  GG1, HH1,  m, p2, q2);
%                 [B_1, B_2, Cm_2_minus_alph, Dm_2_minus_alph]  = coefficientsOfContinuousSpectrum_forPsi2(typeOfCylinder, k_0, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, p2, q2, psi_alph,  1);
                        
        end
        
        %if(k==i && i==2)
        %P_plus_sum = P_plus_sum + a_plus_alph.*conj(a_plus_bet).*...
           %(Cm_2_plus_alph.*conj(Cm_2_plus_bet) + Dm_2_plus_alph.*conj(Dm_2_plus_bet));
%         P_plus_sum = P_plus_sum + Cm_2_plus_alph.*conj(Cm_2_plus_bet) + Dm_2_plus_alph.*conj(Dm_2_plus_bet);
       
        Phase_exp_cs = exp(1i * k_0 * abs(p) * abs(L));
        P_minus_sum = P_minus_sum + (a_minus_alph.*Phase_exp_cs + b_minus_alph).*conj(a_minus_bet.*Phase_exp_cs + b_minus_bet).*...
            (Cm_2_minus_alph.*conj(Cm_2_minus_bet) + Dm_minus_2_alph.*conj(Dm_minus_2_bet)).*psi_alph.*conj(psi_bet);
%         P_minus_sum = P_minus_sum + (Cm_2_minus_alph.*conj(Cm_2_minus_bet) + Dm_minus_2_alph.*conj(Dm_minus_2_bet)).*psi_alph.*conj(psi_bet);
       % end
        
    end
end

% D_plus  = c/(8*pi*k_0^2).*(cot(teta1)).^2.*(c_cs1.*conj(c_cs1) + c_cs2.*conj(c_cs2));
D_plus  = c/(8*pi*k_0^2).*(p./q).^2.*(c_cs1.*conj(c_cs1) + c_cs2.*conj(c_cs2));% (пока только для
%положительного направления оси z) %
% D_minus = c/(2*pi*k_0^2).*(cot(teta2)).^2.*P_minus_sum; % для отрицательного направления оси z
D_minus = c/(2*pi*k_0^2).*(p./q).^2.*P_minus_sum;
D_max = 1;%max(max(abs(D_plus)),max(abs(D_minus)));

% figure(1)
% %subplot(1,2,1)
% plot(teta1/pi,real(D_plus/D_max), teta1/pi, imag(D_plus/D_max),'--');  hold on;
% xlabel('teta (\pi rad)'); ylabel('S_r/max(|S_r|)');
% %subplot(1,2,2)
% plot(teta2/pi,real(D_minus/D_max), teta2/pi, imag(D_minus/D_max),'--'); hold on;
% %xlabel('teta')

% teta = [teta1  teta2];
% D1 = [D_plus zeros(1,size(teta2,2))];
% D2 = [zeros(1,size(teta1,2)) D_minus];
% D = (D1 + D2)/D_max;
% figure(1)
% plot(teta/pi, real(D), 'r'); hold on;

figure(8)
% plot(teta1/pi, D_plus./D_max, 'b'); hold on
plot(q, D_plus./D_max, 'b'); hold on
figure(9)
% plot(teta2/pi, D_minus./D_max, 'r'); hold on
plot(q, D_minus./D_max, 'r'); hold on
xlabel('teta (\pi rad)'); ylabel('S_r/max(|S_r|)');

% Построение в полярных координатах
teta1 = asin(q);
teta2 = pi - asin(q(end:-1:1));

teta = [teta1  teta2];
D1 = [D_plus zeros(1,size(teta2,2))];
D2 = [zeros(1,size(teta1,2)) D_minus(end:-1:1)];
D = (D1 + D2)/D_max;
figure(2)
polar(teta1,(D_plus),'r'); hold on
polar(teta2,(D_minus(end:-1:1)),'r');
% polar(teta, D,'k'); hold on
% polar(-teta, D,'k');









