clear all
clc
tic
systemParameters
toOhms = 898755178736.818;

z0 = 0;
q_0 = sqrt(1-p_0^2);
q_0 = q_0.* (2*(imag(q_0) <= 0)-1);

typeOfCylinder = 'Gyrotropic';
% typeOfCylinder = 'Isotropic';

rho = [0.12:0.01:2] * a_0;
% %%%%%%%%%%%% test field
% Ez = field_of_discreteMode('Ez', typeOfCylinder, rho, q_0, p_0, 1, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, 0);
% figure(5)
% hold on
% plot(rho/a_0, real(Ez),['b' 'o-'],'LineWidth',1)
% plot(rho/a_0, imag(Ez),['r' 'o-'] ,'LineWidth',1)
% hold off
% title('Ez')
% box on
% %%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%% propagation constants of eigenwaves %%%%%%%%%%%%%%%%%%%%%%%%
if(strcmp(typeOfCylinder, 'Isotropic'))
    p_n__of_descreteMode_of_isotropicCyl
elseif(strcmp(typeOfCylinder, 'Gyrotropic'))
    p_n__of_descreteMode_of_gyrotropicCyl
end

% p_n = p_n(real(p_n) < 300);

q_n = sqrt(1-p_n.^2);
q_n = q_n.* (2*(imag(q_n) <= 0)-1);


%%%%%%%%%%%%%%% the calculation of excitation coefficients of the discrete
%%%%%%%%%%%%%%% modes and resistance of waves
for in = 1:size(q_n,1)
    a_smn_forw(in,1) = a_smn_of_discreteMode(typeOfCylinder, q_n(in), p_0,  p_n(in), k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, j_f, j_z, d);
    a_smn_back(in,1) = a_smn_of_discreteMode(typeOfCylinder, q_n(in), p_0, -p_n(in), k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, j_f, j_z, d);
    
    R_of_discreteMode(in,1) =(P_of_discreteMode(typeOfCylinder, q_n(in), p_n(in), k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, a_smn_forw(in,1))+...
                             -P_of_discreteMode(typeOfCylinder, q_n(in),-p_n(in), k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, a_smn_back(in,1))).*...
                              2 / I_0^2 * toOhms;
                          
    Ez(in,1)   = field_of_discreteMode('Ez',   typeOfCylinder, a_0, q_n(in), -p_n(in), 1, k_0, a_0, EE1, -GG1, HH1, MU1, EE, MU, -m, 0);
    Ephi(in,1) = field_of_discreteMode('Ephi', typeOfCylinder, a_0, q_n(in), -p_n(in), 1, k_0, a_0, EE1, -GG1, HH1, MU1, EE, MU, -m, 0);
end
figure(1)
hold on
bar(p_n, R_of_discreteMode)
hold off

Rn_sum = sum(R_of_discreteMode)

 figure
 hold on 
 plot(p_n, abs(Ez./Ephi))
 hold off

% [R_alp1, R_alp2, R_E_freespace, R_H_freespace] = sourceResistance_in_presence_waveguide_and_freeSpace(w_0, typeOfCylinder,...
%         H0, w_H, typeOfmedia, 0, 0, p_0, a_0, MU1, EE, MU, c, m, j_f, j_z, I_0, d, toOhms)

%%
%%%%%%%%%%%%%%% the calculation of dependece of Rsum on d of source

d = [1:1:20] * a_0;

delta = 0;
pp = @(q) sqrt(1 - q.^2);
N_p1 = @(q, p) (-1)^(m+1) * c * p./ (k_0.^2 * q);
N_p2 = @(q, p) (-1)^(m+2) * c * p./ (k_0.^2 * q);

for id=1:size(d,2)
    j_f = 1e7 / (2 * d(id));
    j_z = 0;
    I_0 = sqrt((j_z * 2 * pi * a_0)^2 + (j_f * 2 * d(id))^2);
    
    for in = 1:size(q_n,1)
        a_smn_forw(in,1) = a_smn_of_discreteMode(typeOfCylinder, q_n(in), p_0,  p_n(in), k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, j_f, j_z, d(id));
        a_smn_back(in,1) = a_smn_of_discreteMode(typeOfCylinder, q_n(in), p_0, -p_n(in), k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, j_f, j_z, d(id));
        
        R_of_discreteMode(in,1) =(P_of_discreteMode(typeOfCylinder, q_n(in), p_n(in), k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, a_smn_forw(in,1))+...
            -P_of_discreteMode(typeOfCylinder, q_n(in),-p_n(in), k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, a_smn_back(in,1))).*...
            2 / I_0^2 * toOhms;
    end
    
    Rn_sum_for_j_f(id) = sum(R_of_discreteMode);
    
    %%%%% for case of free space
    P_Ewave_plus  =  -1/4 * quadgk(@(q) abs(a_sma1_freeSpace(q, pp(q), k_0, p_0, a_0, m, c, j_z, j_f, d(id))).^2.*...
        N_p1(q, pp(q)), delta, 1);
    P_Ewave_minus  =   1/4 * quadgk(@(q) abs(a_sma1_freeSpace(q,-pp(q), k_0, p_0, a_0, m, c, j_z, j_f, d(id))).^2.*...
        N_p1(q,-pp(q)), delta, 1);
    
    P_Hwave_plus  =  1/4 * quadgk(@(q) abs(a_sma2_freeSpace(q, pp(q), k_0, p_0, a_0, m, c, j_z, j_f, d(id))).^2.*...
        N_p2(q, pp(q)), delta, 1);
    P_Hwave_minus  =  -1/4 * quadgk(@(q) abs(a_sma2_freeSpace(q,-pp(q), k_0, p_0, a_0, m, c, j_z, j_f, d(id))).^2.*...
        N_p2(q,-pp(q)), delta, 1);
    
%     %%%%% for waveguide
%     P_alp1_forw  =  quadgk(@(q) P_of_continuousWaves_alpha1(typeOfCylinder,...
%         q, p_0, pp(q), k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, 0, j_f, j_z, d(id)), delta, 1);
%     P_alp1_back  = -quadgk(@(q) P_of_continuousWaves_alpha1(typeOfCylinder,...
%         q, p_0,-pp(q), k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, 0, j_f, j_z, d(id)), delta, 1);
%     
%     P_alp2_forw  = -quadgk(@(q) P_of_continuousWaves_alpha2(typeOfCylinder,...
%         q, p_0, pp(q), k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, 0, j_f, j_z, d(id)), delta, 1);
%     P_alp2_back  = quadgk(@(q) P_of_continuousWaves_alpha2(typeOfCylinder,...
%         q, p_0,-pp(q), k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, 0, j_f, j_z, d(id)), delta, 1);
    
    R_E_freespace_f(id) = abs((P_Ewave_plus + P_Ewave_minus)) * 2 / I_0^2 * toOhms;
%     R_alp1        = abs((P_alp1_forw + P_alp1_back))  * 2 / I_0^2   * toOhms
    R_H_freespace_f(id) = abs((P_Hwave_plus + P_Hwave_minus)) * 2 / I_0^2 * toOhms;
%     R_alp2        = abs((P_alp2_forw + P_alp2_back))  * 2 / I_0^2   * toOhms
    
    
    
    
    j_f = 0;
    j_z = 2e6 / (2 * pi * a_0);
    I_0 = sqrt((j_z * 2 * pi * a_0)^2 + (j_f * 2 * d(id))^2);    
    for in = 1:size(q_n,1)
        a_smn_forw(in,1) = a_smn_of_discreteMode(typeOfCylinder, q_n(in), p_0,  p_n(in), k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, j_f, j_z, d(id));
        a_smn_back(in,1) = a_smn_of_discreteMode(typeOfCylinder, q_n(in), p_0, -p_n(in), k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, j_f, j_z, d(id));
        
        R_of_discreteMode(in,1) =(P_of_discreteMode(typeOfCylinder, q_n(in), p_n(in), k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, a_smn_forw(in,1))+...
            -P_of_discreteMode(typeOfCylinder, q_n(in),-p_n(in), k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, a_smn_back(in,1))).*...
            2 / I_0^2 * toOhms;
    end
    
    Rn_sum_for_j_z(id) = sum(R_of_discreteMode);
    
       
    %%%%% for case of free space
    P_Ewave_plus  =  -1/4 * quadgk(@(q) abs(a_sma1_freeSpace(q, pp(q), k_0, p_0, a_0, m, c, j_z, j_f, d(id))).^2.*...
        N_p1(q, pp(q)), delta, 1);
    P_Ewave_minus  =   1/4 * quadgk(@(q) abs(a_sma1_freeSpace(q,-pp(q), k_0, p_0, a_0, m, c, j_z, j_f, d(id))).^2.*...
        N_p1(q,-pp(q)), delta, 1);
    
    P_Hwave_plus  =  1/4 * quadgk(@(q) abs(a_sma2_freeSpace(q, pp(q), k_0, p_0, a_0, m, c, j_z, j_f, d(id))).^2.*...
        N_p2(q, pp(q)), delta, 1);
    P_Hwave_minus  =  -1/4 * quadgk(@(q) abs(a_sma2_freeSpace(q,-pp(q), k_0, p_0, a_0, m, c, j_z, j_f, d(id))).^2.*...
        N_p2(q,-pp(q)), delta, 1);
    
%     %%%%% for waveguide
%     P_alp1_forw  =  quadgk(@(q) P_of_continuousWaves_alpha1(typeOfCylinder,...
%         q, p_0, pp(q), k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, 0, j_f, j_z, d(id)), delta, 1);
%     P_alp1_back  = -quadgk(@(q) P_of_continuousWaves_alpha1(typeOfCylinder,...
%         q, p_0,-pp(q), k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, 0, j_f, j_z, d(id)), delta, 1);
%     
%     P_alp2_forw  = -quadgk(@(q) P_of_continuousWaves_alpha2(typeOfCylinder,...
%         q, p_0, pp(q), k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, 0, j_f, j_z, d(id)), delta, 1);
%     P_alp2_back  = quadgk(@(q) P_of_continuousWaves_alpha2(typeOfCylinder,...
%         q, p_0,-pp(q), k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, 0, j_f, j_z, d(id)), delta, 1);
    
    R_E_freespace_z(id) = abs((P_Ewave_plus + P_Ewave_minus)) * 2 / I_0^2 * toOhms;
%     R_alp1        = abs((P_alp1_forw + P_alp1_back))  * 2 / I_0^2   * toOhms
    R_H_freespace_z(id) = abs((P_Hwave_plus + P_Hwave_minus)) * 2 / I_0^2 * toOhms;
%     R_alp2        = abs((P_alp2_forw + P_alp2_back))  * 2 / I_0^2   * toOhms
end
figure(45)
hold on
plot(d/a_0, Rn_sum_for_j_f,'b')
plot(d/a_0, Rn_sum_for_j_z,'r')
plot(d/a_0, (R_E_freespace_f + R_H_freespace_f),'b*')
plot(d/a_0, (R_E_freespace_z + R_H_freespace_z),'r*')
hold off



%% sum of resistance of modes
w0s_min = w_0;
d = a_0;
DD =  1/2.5;
w0s_max = DD * w_0;
w0s = [w0s_min:0.1*(w0s_max-w0s_min):w0s_max];

%%%%%% resistance of mode
[R_of_discreteMode, p_n_new, Pin, Pcin] = sourceResistance_of_modes(w0s, p_n, typeOfCylinder,...
                H0, w_H, typeOfmedia, 0, 0, p_0, a_0, c, m, j_f, j_z, I_0, d, toOhms);
%%%%%% sum of resistance of mode
 sumR_of_discreteMode = sum(R_of_discreteMode,1);

figure(4)
hold on
bar(p_n, R_of_discreteMode(:,1))
hold off

II = size(p_n_new,2);
% for II=1:size(p_n_new,2)
figure(60+II)
hold on
bar(p_n, R_of_discreteMode(:,II))
hold off
% end

figure(6)
hold on
plot(w0s/w_H,real(Pin),'g--',w0s/w_H,Pcin,'r--')
for J=1:size(p_n_new,1)
      plot(w0s/w_H,p_n_new(J,:))
      text(w0s(1)/w_H, p_n_new(J,1),['n=' num2str(J)]);
end
hold off


figure(7)
hold on
plot(w0s/w_H, sumR_of_discreteMode)
hold off






%%%%%%%%%%%%%%% the calculation of excitation coefficients of the continuous-spectrum waves
% upper_Bound = 1000;
% N = 100;
% [dq_simp, q, p] = quadratureMethod_forIntegralEqs('simpson', N, q_0, 1e-8, upper_Bound);

% p = sqrt(1 - q.^2);
% p = real(p) - 1i * abs(imag(p));
% 
% a_sma1_forw = a_sma_of_continuousWaves_alpha1(typeOfCylinder, q, p_0, p, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, 0, I_f, I_z, d);
% a_sma1_back = a_sma_of_continuousWaves_alpha1(typeOfCylinder, q, p_0,-p, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, 0, I_f, I_z, d);

%% calculation of source resistance in free space and in the presence of the waveguide
delta = 0;
pp = @(q) sqrt(1 - q.^2);
N_p1 = @(q, p) (-1)^(m+1) * c * p./ (k_0.^2 * q);
N_p2 = @(q, p) (-1)^(m+2) * c * p./ (k_0.^2 * q);

%%%%% for free space
P_Ewave_plus  =  -1/4 * quadgk(@(q) abs(a_sma1_freeSpace(q, pp(q), k_0, p_0, a_0, m, c, j_z, j_f, d)).^2.*...
    N_p1(q, pp(q)), delta, 1);
P_Ewave_minus  =   1/4 * quadgk(@(q) abs(a_sma1_freeSpace(q,-pp(q), k_0, p_0, a_0, m, c, j_z, j_f, d)).^2.*...
    N_p1(q,-pp(q)), delta, 1);

P_Hwave_plus  =  1/4 * quadgk(@(q) abs(a_sma2_freeSpace(q, pp(q), k_0, p_0, a_0, m, c, j_z, j_f, d)).^2.*...
    N_p2(q, pp(q)), delta, 1);
P_Hwave_minus  =  -1/4 * quadgk(@(q) abs(a_sma2_freeSpace(q,-pp(q), k_0, p_0, a_0, m, c, j_z, j_f, d)).^2.*...
    N_p2(q,-pp(q)), delta, 1);

%%%%% for waveguide
P_alp1_forw  =  quadgk(@(q) P_of_continuousWaves_alpha1(typeOfCylinder,...
    q, p_0, pp(q), k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, 0, j_f, j_z, d), delta, 1);
P_alp1_back  = -quadgk(@(q) P_of_continuousWaves_alpha1(typeOfCylinder,...
    q, p_0,-pp(q), k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, 0, j_f, j_z, d), delta, 1);

P_alp2_forw  = -quadgk(@(q) P_of_continuousWaves_alpha2(typeOfCylinder,...
    q, p_0, pp(q), k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, 0, j_f, j_z, d), delta, 1);
P_alp2_back  = quadgk(@(q) P_of_continuousWaves_alpha2(typeOfCylinder,...
    q, p_0,-pp(q), k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, 0, j_f, j_z, d), delta, 1);

R_E_freespace = abs((P_Ewave_plus + P_Ewave_minus)) * 2 / I_0^2 * toOhms
R_alp1        = abs((P_alp1_forw + P_alp1_back))  * 2 / I_0^2   * toOhms
R_H_freespace = abs((P_Hwave_plus + P_Hwave_minus)) * 2 / I_0^2 * toOhms
R_alp2        = abs((P_alp2_forw + P_alp2_back))  * 2 / I_0^2   * toOhms



% a_sma2_forw = a_sma_of_continuousWaves_alpha2(typeOfCylinder, q, p_0, p, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, 0, I_f, I_z, d);
% a_sma2_back = a_sma_of_continuousWaves_alpha2(typeOfCylinder, q, p_0,-p, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, 0, I_f, I_z, d);
% 
% figure(3)
% hold on
% plot(q, real(a_sma1_forw), q, real(a_sma1_back), 'r')
% hold off

% figure(4)
% hold on
% plot(q, real(a_sma2_forw), q, real(a_sma2_back), 'r')
% hold off





%% calculation of dependece of source resistance in free space and in the presence of the waveguide on frequency

% w0s_min = w_0;
% DD =  1/2.5;
% w0s_max = DD * w_0;
% w0s = [w0s_min:0.1*(w0s_max-w0s_min):w0s_max];

[R_alp1, R_alp2, R_E_freespace, R_H_freespace] = sourceResistance_in_presence_waveguide_and_freeSpace(w0s, typeOfCylinder,...
        H0, w_H, typeOfmedia, 0, 0, p_0, a_0, MU1, EE, MU, c, m, j_f, j_z, I_0, d, toOhms);

figure(2)
   hold on
   plot(w0s/w_H,R_E_freespace,'b',w0s/w_H,R_H_freespace,'r')
   plot(w0s/w_H,R_alp1,'b--',w0s/w_H,R_alp2,'r--')
   hold off
   
   NofmaxB = 1;
   text(w0s(NofmaxB)/w_H, abs(R_E_freespace(NofmaxB)),'R_{Efreespace}');
   text(w0s(NofmaxB)/w_H, abs(R_H_freespace(NofmaxB)),'R_{Hfreespace}');
   text(w0s(NofmaxB)/w_H, abs(R_alp1(NofmaxB)),'R_{alp1}');
   text(w0s(NofmaxB)/w_H, abs(R_alp2(NofmaxB)),'R_{alp2}');
hold off

%% dependences of Rn_sum and Rcont on d
tic
clear R_of_discreteMode

I = 1;
dd = [d:0.1 * d:10*d];
for idd = dd
    for in = 1:size(q_n,1)
        a_smn_forw = a_smn_of_discreteMode(typeOfCylinder, q_n(in), p_0,  p_n(in), k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, j_f, j_z, idd);
        a_smn_back = a_smn_of_discreteMode(typeOfCylinder, q_n(in), p_0, -p_n(in), k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, j_f, j_z, idd);
        
        R_of_discreteMode(in,1) =(P_of_discreteMode(typeOfCylinder, q_n(in), p_n(in), k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, a_smn_forw)+...
            -P_of_discreteMode(typeOfCylinder, q_n(in),-p_n(in), k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, a_smn_back)).*...
            2 / I_0^2 * toOhms;
    end
Rn_sum(I) = sum(R_of_discreteMode);

[R_alp1(I), R_alp2(I), R_E_freespace, R_H_freespace] = sourceResistance_in_presence_waveguide_and_freeSpace(w_0, typeOfCylinder,...
        H0, w_H, typeOfmedia, 0, 0, p_0, a_0, MU1, EE, MU, c, m, j_f, j_z, I_0, idd, toOhms);

I = I + 1;
end


figure(1)
hold on
  plot(dd, Rn_sum)
  plot(dd, R_alp1, 'b--')
  plot(dd, R_alp2, 'r--')
hold off


toc
  
  
