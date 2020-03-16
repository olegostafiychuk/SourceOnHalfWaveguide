clearvars
clc
tic
systemParameters
toOhms = 898755178736.818;

z0 = 0;
% q_0 = sqrt(1-p_0^2);
% q_0 = q_0.* (2*(imag(q_0) <= 0)-1);

typeOfCylinder = 'Gyrotropic';
% typeOfCylinder = 'Isotropic';

%%%%%%%%%%%%%% propagation constants of eigenwaves %%%%%%%%%%%%%%%%%%%%%%%%
if(strcmp(typeOfCylinder, 'Isotropic'))
    p_n__of_descreteMode_of_isotropicCyl
elseif(strcmp(typeOfCylinder, 'Gyrotropic'))
%     p_n__of_descreteMode_of_gyrotropicCyl
    p_n__of_descreteMode_of_gyrotropicCyl
end
% p_n = p_n(real(p_n) < 1000);
p_n = p_n(1:100);

q_n = sqrt(1-p_n.^2);
q_n = q_n.* (2*(imag(q_n) <= 0)-1);

% EE1 = 1.0000001;
% GG1 = 0.0000001;
% HH1 = 1.0000001;

%control of boundary conditions and ratios of the components
% rho = [0:0.0001:1.05] * a_0;
% %%%%%%%%%%%% test field
% Ez = field_of_discreteMode('Ephi', typeOfCylinder, rho, q_n(8), p_n(8), 1, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, 0);
% 
% figure(5)
% hold on
% plot(rho/a_0, real(Ez),'b','LineWidth',1)
% plot(rho/a_0, imag(Ez),'r','LineWidth',1)
% hold off
% title('Ez')
% box on


% Ez_v = field_of_discreteMode('Ez', typeOfCylinder, a_0, q_n, p_n, 1, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, 0);
% Ephi_v = field_of_discreteMode('Ephi', typeOfCylinder, a_0, q_n, p_n, 1, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, 0);
% bar(p_n, abs(Ez_v)./abs(Ephi_v));

%% %%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%% the calculation of excitation coefficients of the discrete
%%%%%%%%%%%%%%% modes and resistance for the discrete-spectrum waves

j_z_vec = 0.254*j_f;
% j_z = 2e6 / (2 * pi * a_0);
% j_z_vec = linspace(0,10).*j_f;

p_0 = 0;
% p_0_vec = 0:1:200;

% % Nu_e = 0;
% clear p_n; clear q_n;
% p_n_losses; %этот скрипт подт€гивает значени€ продольных волновых чисел и столбец столкновений из файла disperson_losses_m1.mat, наход€щегос€ в этой же директории
% p_n_los = p_n_real + 1i*p_n_imag;
% p_n = p_n(real(p_n) < 1000);
% Nu_e_vec = Nu_e_vec(71); %<-- здесь задаем число от 1 до MAX - выбираем необходимое знаечение столкновений (по возрастанию)

I = 0;
for j_z = j_z_vec
% for p_0 = p_0_vec
% for Nu_e = Nu_e_vec
    I = I + 1
       
%     p_n = p_n_los(:,I); %<-- здесь задаем число от 1 до 51 - это будет столбец p_n дл€ соотвтетствующего значени€ столкновений (дл€ всех столбцов ставим I)
%     q_n = sqrt(1-p_n.^2);
%     q_n = q_n.* (2*(imag(q_n) <= 0)-1);
   
    [EE1, GG1, HH1, c] = channelparameters_sources(H0, typeOfmedia, w_0, 0, 0, n_e, Nu_e);
            waveguideParameters.EEinner = EE1;
            waveguideParameters.GGinner = GG1;
            waveguideParameters.HHinner = HH1;
            
    I_0 = sqrt((j_z * 2 * pi * a_0)^2 + (j_f * 2 * d)^2);  
    
for in = 1:size(q_n,1)   
             
        a_smn_forw(in) = a_smn_of_discreteMode(typeOfCylinder, q_n(in), p_0,  p_n(in), k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, j_f, j_z, d);
        a_smn_back(in) = a_smn_of_discreteMode(typeOfCylinder, q_n(in), p_0, - p_n(in), k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, j_f, j_z, d);
        N_n(in) = Norm_of_descreteMode(typeOfCylinder, q_n(in), p_n(in), k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m);
        
        %function for calculating of power modified by Vasiliy
        P_n_plus(in) = P_of_discreteMode_NEW( p_n(in), waveguideParameters, sourceParameters, a_smn_forw(in));
        P_n_minus(in) = P_of_discreteMode_NEW( - p_n(in), waveguideParameters, sourceParameters, a_smn_back(in));
        
        %Power_int_plus(in) = powerOfDiscreteSpectrumMode(p_n(in), waveguideParameters, sourceParameters);
        
        P_n(in) = P_n_plus(in) - P_n_minus(in);
        R_n(in) = P_n(in).*2 / I_0.^2 * toOhms;
    

end

%%%%% for continuous-spectrum waves
    clear P_alp1_forw P_alp1_back P_alp2_forw P_alp2_back
    delta = 0;
    pp = @(q) sqrt(1 - q.^2);
    P_alp1_forw  =  quadgk(@(q) P_of_continuousWaves_alpha1(typeOfCylinder,...
        q, p_0, pp(q), k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, 0, j_f, j_z, d), delta, 1);
    P_alp1_back  = quadgk(@(q) P_of_continuousWaves_alpha1(typeOfCylinder,...
        q, p_0,- pp(q), k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, 0, j_f, j_z, d), delta, 1);

    P_alp2_forw  = quadgk(@(q) P_of_continuousWaves_alpha2(typeOfCylinder,...
        q, p_0, pp(q), k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, 0, j_f, j_z, d), delta, 1);
    P_alp2_back  = quadgk(@(q) P_of_continuousWaves_alpha2(typeOfCylinder,...
        q, p_0,- pp(q), k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, 0, j_f, j_z, d), delta, 1);

    P_cs(I) = abs(P_alp1_forw + P_alp2_forw) + abs(P_alp1_back + P_alp2_back); 
    R_cs(I) = P_cs(I).*2 / I_0.^2 * toOhms;
    
    P_mod(I) = sum(P_n);
    R_mod(I) = sum(R_n);

%     R_cs(ie) = P_cs(ij).*2 / I_0.^2 * toOhms;
    
%     figure(11)
%     bar(real(p_n), P_n.*2 / I_0.^2 * toOhms);
%     figure(12)
%     bar(real(p_n), abs(a_smn_forw).^2);
%     figure(13)
%     bar(real(p_n), Power_int_plus);
%     figure(14)
%     bar(real(p_n), N_n./Power_int_plus);
%     figure(15)
%     bar(real(p_n), imag(N_n)); title('Im N_{n}')
end
figure(13)
plot(j_z_vec/j_f, R_mod,'k--'); hold on;
figure(14)
plot(j_z_vec/j_f, R_cs,'k--'); hold on;
% plot(p_0_vec, R_mod,'k-.'); hold on
% figure(16)
% plot(Nu_e_vec/w_H, R_mod,'k-*'); hold on
% xlabel('\nu_{e}/\omega_H');
% ylabel('R_{mod}');
% figure(17)
% plot(Nu_e_vec/w_H, R_cs,'k-*'); hold on
% xlabel('\nu_{e}/\omega_H');
% ylabel('R_{cs}');
toc

%% calculation of source resistance in free space
delta = 0;
pp = @(q) sqrt(1 - q.^2);
N_p1 = @(q, p) (-1)^(m+1) * c * p./ (k_0.^2 * q);
N_p2 = @(q, p) (-1)^(m+2) * c * p./ (k_0.^2 * q);

q=1e-2:1e-3:1;
p = sqrt(1 - q.^2);

P_plus_alph1 = P_of_continuousWaves_alpha1(typeOfCylinder, q, p_0, p, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, 0, j_f, j_z, d);
P_minus_alph1 = - P_of_continuousWaves_alpha1(typeOfCylinder, q, p_0, - p, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, 0, j_f, j_z, d);
R_alp1_q = abs(P_plus_alph1 + P_minus_alph1)*2 / I_0^2 * toOhms;

P_plus_alph2 = P_of_continuousWaves_alpha2(typeOfCylinder, q, p_0, p, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, 0, j_f, j_z, d);
P_minus_alph2 = - P_of_continuousWaves_alpha2(typeOfCylinder, q, p_0, - p, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, 0, j_f, j_z, d);
R_alp2_q = abs(P_plus_alph2 + P_minus_alph2)*2 / I_0^2 * toOhms;

figure(60)
plot(q,R_alp1_q,'r',q,R_alp2_q,'b');
%set(gca,'YScale','log');

%% %%% for free space
P_Ewave_plus  =  -1/4 * quadgk(@(q) abs(a_sma1_freeSpace(q, pp(q), k_0, p_0, a_0, m, c, j_z, j_f, d)).^2.*...
    N_p1(q, pp(q)), delta, 1);
P_Ewave_minus  =   1/4 * quadgk(@(q) abs(a_sma1_freeSpace(q,-pp(q), k_0, p_0, a_0, m, c, j_z, j_f, d)).^2.*...
    N_p1(q,-pp(q)), delta, 1);

P_Hwave_plus  =  1/4 * quadgk(@(q) abs(a_sma2_freeSpace(q, pp(q), k_0, p_0, a_0, m, c, j_z, j_f, d)).^2.*...
    N_p2(q, pp(q)), delta, 1);
P_Hwave_minus  =  -1/4 * quadgk(@(q) abs(a_sma2_freeSpace(q,-pp(q), k_0, p_0, a_0, m, c, j_z, j_f, d)).^2.*...
    N_p2(q,-pp(q)), delta, 1);

R_E_freespace = abs((P_Ewave_plus + P_Ewave_minus)) * 2 / I_0^2 * toOhms;
R_H_freespace = abs((P_Hwave_plus + P_Hwave_minus)) * 2 / I_0^2 * toOhms;
R_freespace = R_E_freespace + R_H_freespace

%%% рамка с током
R_ring = 2*pi^2/(3*c)*(k_0*a_0)^4.*toOhms
R_dip = 2/(3*c)*(k_0*2*d)^2.*toOhms
  
  
