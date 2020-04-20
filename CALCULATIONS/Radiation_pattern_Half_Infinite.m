clearvars
clc

% вычисление диаграммы направленности несимметричного источника
typeOfCylinder = 'Gyrotropic';
% typeOfCylinder = 'Isotropic';
systemParameters

saveVarsMat = load('C:\Users\ќлег\Dropbox\PAPERS\URSI GASS 2020\coefficients_0025wH_n1e13_L_1_5lambda0.mat');
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
q = q(1:ind);
p = p(1:ind);
a_cs1_minus = a_cs1_minus(1:ind);
a_cs2_minus = a_cs2_minus(1:ind);
b_cs1_minus = b_cs1_minus(1:ind);
b_cs2_minus = b_cs2_minus(1:ind);
c_cs1 = c_cs1(1:ind);
c_cs2 = c_cs2(1:ind);


% plot(q,real(a_cs2_minus),'r-x',q,imag(a_cs2_minus),'b-o');

% EE1 =  1.00001;
% GG1 = 0.000001;
% HH1 =  1.00001;

% step = 1/100;
% teta1 = (0.0001:step:0.49999)*pi;
% teta2 = (0.50001:step:0.9999)*pi;
% %teta2 = teta1;
% q1 = sin(teta1);
% q2 = sin(teta2);
% p1 = sqrt(1 - q1.^2);%cos(teta1);%
% p2 = sqrt(1 - q2.^2);%cos(teta2);%
% z = 0;
% psi1   = psi1_q__gyrotropic(k_0, k_0, a_0, EE1,  GG1, HH1,  m, p2, q2);
% psi2   = psi2_q__gyrotropic(k_0, k_0, a_0, EE1,  GG1, HH1,  m, p2, q2);

%%%%%%%%%%%%%%%%%%%%%
% if(strcmp(typeOfCylinder, 'Isotropic'))
%     p_n__of_descreteMode_of_isotropicCyl
% elseif(strcmp(typeOfCylinder, 'Gyrotropic'))
%     p_n__of_descreteMode_of_gyrotropicCyl
% %     p_n = p_n(1:31);
%     p_n(7:end) = - p_n(7:end);
% end
% 
% q_n = sqrt(1-p_n.^2);
% q_n = q_n.* (2*(imag(q_n) <= 0)-1);
% 
% for in = 1:size(q_n,1)
%     a_smn_forw(in,1) = a_smn_of_discreteMode_Decorator(typeOfCylinder, q_n(in),  p_n(in), waveguideParameters, sourceParameters);
%     a_smn_back(in,1) = a_smn_of_discreteMode_Decorator(typeOfCylinder, q_n(in), -p_n(in), waveguideParameters, sourceParameters);
% end
% 
% for iq = 1:size(q_cs,2)
%     a_cs_alfa1_forw(iq) = a_sma_of_continuousWaves_alpha1_Decorator(typeOfCylinder, q_cs(iq),  p_cs(iq), waveguideParameters, sourceParameters);
%     a_cs_alfa1_back(iq) = a_sma_of_continuousWaves_alpha1_Decorator(typeOfCylinder, q_cs(iq), -p_cs(iq), waveguideParameters, sourceParameters);
%     a_cs_alfa2_forw(iq) = a_sma_of_continuousWaves_alpha2_Decorator(typeOfCylinder, q_cs(iq),  p_cs(iq), waveguideParameters, sourceParameters);
%     a_cs_alfa2_back(iq) = a_sma_of_continuousWaves_alpha2_Decorator(typeOfCylinder, q_cs(iq), -p_cs(iq), waveguideParameters, sourceParameters);
% end
% 
% LL = 4*d;%lambda_0;%[0.2 0.5 1.5]*lambda_0;
% L = - LL(il);
% sourceParameters.zCoordinate = - LL(il);
% 
% continuousSperctrumCharacters.dq_simp = dq_simp;
% continuousSperctrumCharacters.q = q_cs;
% continuousSperctrumCharacters.p = p_cs;
% continuousSperctrumCharacters.a_cs_alfa1_forw = a_cs_alfa1_forw;
% continuousSperctrumCharacters.a_cs_alfa2_forw = a_cs_alfa2_forw;
% 
% discreteSperctrumCharacters.q_n = q_n;
% discreteSperctrumCharacters.p_n = p_n;
% discreteSperctrumCharacters.a_smn_forw = a_smn_forw;
% 
% [b_p_field_1_back, b_p_field_2_back, b_n_back, a_p_Ebeam_forw,a_p_Hbeam_forw] =...           
%           scatteringCoeffs_of_EigenWavesOfHalfInfinityGyrotropCylEdge(typeOfCylinder,...
%                continuousSperctrumCharacters, discreteSperctrumCharacters,waveguideParameters, sourceParameters);

%%%%%%%%%%%%%%%%%%%%

teta1 = asin(q);
teta2 = pi/2 + asin(q);
q1 = sin(teta1);
q2 = sin(teta2);
p1 = cos(teta1);%sqrt(1 - q1.^2);%cos(teta1);
p2 = cos(teta2);%sqrt(1 - q2.^2);%
% plot(teta1, p1); hold on;
% plot(teta2, p2)
a_cs1_minus = a_cs1_minus(end:-1:1);
a_cs2_minus = a_cs2_minus(end:-1:1);
b_cs1_minus = b_cs1_minus(end:-1:1);
b_cs2_minus = b_cs2_minus(end:-1:1);

P_plus_sum = zeros(1,size(teta1,2));
P_minus_sum = zeros(1,size(teta2,2));
% P_plus_sum = zeros(1,size(q1,2));
% P_minus_sum = zeros(1,size(q2,2));
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
                psi_bet = psi1_q__gyrotropic(k_0, k_0, a_0, EE1,  GG1, HH1,  m, p2, q2);
                [B_1, B_2, Cm_2_minus_bet, Dm_minus_2_bet]  = coefficientsOfContinuousSpectrum(typeOfCylinder, k_0, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, p2, q2, psi_bet,  1);
                                
              
            case 2
                %a_plus_bet = a_sma_of_continuousWaves_alpha2(typeOfCylinder, q1, p_0, p1, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, z, j_f, j_z, d);
%                 a_minus_bet = a_sma_of_continuousWaves_alpha2(typeOfCylinder, q2, p_0, -p2, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, z, j_f, j_z, d);
                %[B_1, B_2, Cm_2_plus_bet, Dm_2_plus_bet]  = coefficientsOfContinuousSpectrum_forPsi2(typeOfCylinder, k_0, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, p1, q1, psi2,  1);
%                 psi_bet = psi2_q__gyrotropic(k_0, k_0, a_0, EE1,  GG1, HH1,  m, -p2, q2);
%                 [B_1, B_2, Cm_2_minus_bet, Dm_minus_2_bet]  = coefficientsOfContinuousSpectrum_forPsi2(typeOfCylinder, k_0, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, -p2, q2, psi_bet,  1);
                
                a_minus_bet = a_cs2_minus;
                b_minus_bet = b_cs2_minus;
                psi_bet = psi2_q__gyrotropic(k_0, k_0, a_0, EE1,  GG1, HH1,  m, p2, q2);
                [B_1, B_2, Cm_2_minus_bet, Dm_minus_2_bet]  = coefficientsOfContinuousSpectrum_forPsi2(typeOfCylinder, k_0, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, p2, q2, psi_bet,  1);
                
                
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
                psi_alph = psi1_q__gyrotropic(k_0, k_0, a_0, EE1,  GG1, HH1,  m, p2, q2);
                [B_1, B_2, Cm_2_minus_alph, Dm_2_minus_alph]  = coefficientsOfContinuousSpectrum(typeOfCylinder, k_0, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, p2, q2, psi_alph,  1);
            
            
            case 2
                %a_plus_alph = a_sma_of_continuousWaves_alpha2(typeOfCylinder, q1, p_0, p1, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, z, j_f, j_z, d);
%                 a_minus_alph = a_sma_of_continuousWaves_alpha2(typeOfCylinder, q2, p_0, -p2, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, z, j_f, j_z, d);
                %[B_1, B_2, Cm_2_plus_alph, Dm_2_plus_alph]  = coefficientsOfContinuousSpectrum_forPsi2(typeOfCylinder, k_0, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, p1, q1, psi2,  1);
%                 psi_alph = psi2_q__gyrotropic(k_0, k_0, a_0, EE1,  GG1, HH1,  m, -p2, q2);
%                 [B_1, B_2, Cm_2_minus_alph, Dm_2_minus_alph]  = coefficientsOfContinuousSpectrum_forPsi2(typeOfCylinder, k_0, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, -p2, q2, psi_alph,  1);
%                 
                a_minus_alph = a_cs2_minus;
                b_minus_alph = b_cs2_minus;
                psi_alph = psi2_q__gyrotropic(k_0, k_0, a_0, EE1,  GG1, HH1,  m, p2, q2);
                [B_1, B_2, Cm_2_minus_alph, Dm_2_minus_alph]  = coefficientsOfContinuousSpectrum_forPsi2(typeOfCylinder, k_0, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, p2, q2, psi_alph,  1);
                        
        end
        
        %if(k==i && i==2)
        %P_plus_sum = P_plus_sum + a_plus_alph.*conj(a_plus_bet).*...
           %(Cm_2_plus_alph.*conj(Cm_2_plus_bet) + Dm_2_plus_alph.*conj(Dm_2_plus_bet));
%         P_plus_sum = P_plus_sum + Cm_2_plus_alph.*conj(Cm_2_plus_bet) + Dm_2_plus_alph.*conj(Dm_2_plus_bet);
       
        Phase_exp_cs = 1;%exp(1i * k_0 * abs(p2) * abs(L));
        P_minus_sum = P_minus_sum + (a_minus_alph.*Phase_exp_cs + b_minus_alph).*conj(a_minus_bet.*Phase_exp_cs + b_minus_bet).*...
            (Cm_2_minus_alph.*conj(Cm_2_minus_bet) + Dm_2_minus_alph.*conj(Dm_minus_2_bet));%*psi_alph.*conj(psi_bet);
%         P_minus_sum = P_minus_sum + Cm_2_minus_alph.*conj(Cm_2_minus_bet) + Dm_2_minus_alph.*conj(Dm_minus_2_bet);
       % end
        
    end
end

D_plus  = c/(8*pi*k_0^2).*(cot(teta1)).^2.*(c_cs1.*conj(c_cs1) + c_cs2.*conj(c_cs2));
% D_plus = c/(2*pi*k_0^2).*(cot(teta1)).^2.*P_plus_sum; % (пока только дл€
%положительного направлени€ оси z) %
D_minus = c/(2*pi*k_0^2).*(cot(teta2)).^2.*P_minus_sum; % дл€ отрицательного направлени€ оси z
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
plot(teta1/pi, D_plus./D_max, 'b'); hold on
% teta_inv = pi.*ones(1,length(teta2)) - teta2;
figure(9)
plot(teta2/pi, D_minus./D_max, 'r'); hold on

xlabel('teta (\pi rad)'); ylabel('S_r/max(|S_r|)');

% ѕостроение в пол€рных координатах
teta = [teta1  teta2];
D1 = [D_plus zeros(1,size(teta2,2))];
D2 = [zeros(1,size(teta1,2)) D_minus];
D = (D1 + D2)/D_max;
figure(2)
polar(teta, D,'k'); hold on
polar(-teta, D,'k');

% polar(teta1,real(D_plus/D_max));  %hold on;
% polar(teta2,real(D_minus/D_max));









