%%%%%% diffraction of bessel beams E and H polarized on edge of cylinder

clear all
tic
systemParameters
% EE1 = 1.0002;

%%%%%% paramters of beam
p_0 = [
%     3.37787970882855;
    30.9884595883368;
          ];
       
q_0 = sqrt(1-p_0^2);
q_0 = q_0.* (2*(imag(q_0) <= 0)-1);

upper_Bound = 600000;%%%%%%%  for experiment of Markov and Gindendurg: n_e = 1e9; H_0 = 0.5; Nu_e = 0; a_0 = 10 cm; w_0 = 1e6;
N = 1024;
[dq_simp, q_back, p] = quadratureMethod_forIntegralEqs('simpson', N, 1i, 1e-8, upper_Bound);
% [dq_simp, q_back, p_cs] = quadratureMethod_forIntegralEqs('difficult_simpsonMarkovGildendurgExpWithCollisions', N, 1i, 1e-7, upper_Bound);



% upper_Bound = 1500000;
% N = 15;
% N_upper = 1000;
% p_n__of_descreteMode_of_gyrotropicCyl;
% q_n00 = (sqrt((1-p_n(abs(imag(p_n))<350000).^2)));
% q_n00 = sqrt((1-p_n.^2));
% q_n0 = [2];
% % for iq0 = 1:size(q_n00,1)
% %     if(q_n00(iq0) == abs(real(q_0)))
% %         q_n0 = [q_n0; q_n00(iq0)-100*abs(imag(q_0)); q_n00(iq0); q_n00(iq0)+100*abs(imag(q_0))];
% %     else
% %         q_n0 = [q_n0; q_n00(iq0)];
% %     end
% % end
% for iq0 = 1:size(q_n00,1)
%     if(abs(imag(q_n00(iq0))) < 1 && abs(imag(q_n00(iq0))) > 0)
%         q_n0 = [q_n0; real(q_n00(iq0))-100*abs(imag(q_n00(iq0))); real(q_n00(iq0)); real(q_n00(iq0))+100*abs(imag(q_n00(iq0)))];
%     else
%         q_n0 = [q_n0; real(q_n00(iq0))];
%     end
% end
% q_n0 = [q_n0; q_n0(end) + 300000];
% [dq_simp, q_back, p] = quadratureMethod_forIntegralEqs_multipleIntervals('difficult_simpsonMarkovGildendurgExpWithCollisions', N, N_upper, q_n0, 1e-8, upper_Bound);


upper_Bound = 1200000;
N = 16;         %%% must be even
N_upper = 1024;  %%% must be even
p_n__of_descreteMode_of_gyrotropicCyl;
% q_n00 = (sqrt((1-p_n(abs(imag(p_n))<150000).^2)));
q_n00 = sqrt((1-p_n.^2));
q_n0 = [2];
q_n0 = [2; 18.0847; 50; 500];
for iq0 = 2:size(q_n00,1)
    if(abs(imag(q_n00(iq0))) < 1 && abs(imag(q_n00(iq0))) > 0)
        q_n0 = [q_n0; real(q_n00(iq0))-100*abs(imag(q_n00(iq0))); real(q_n00(iq0)); real(q_n00(iq0))+100*abs(imag(q_n00(iq0)))];
%         q_n0 = [q_n0; real(q_n00(iq0))];
    else
        q_n0 = [q_n0; real(q_n00(iq0))];
    end
end
q_n0 = [q_n0; q_n0(end) + 300000];
[dq_simp, q_back, p] = quadratureMethod_forIntegralEqs_multipleIntervals('difficult_simpsonMarkovGildendurgExpWithCollisions', N, N_upper, q_n0, 1e-8, upper_Bound);



upper_Bound = 40;
N = 1600;
[dq_simp, q_back, p] = quadratureMethod_forIntegralEqs('simpson', N, 1i, 1e-8, upper_Bound);
% [dq_simp, q_back, p_cs] = quadratureMethod_forIntegralEqs('difficult_simpsonMarkovGildendurgExp', N, 1i, 1e-7, upper_Bound);


%  p_n__of_descreteMode_of_gyrotropicCyl;
 
 
p_back = sqrt(1 - q_back.^2);
p_back = p_back.*(2*(imag(p_back) <= 0)-1);
    
%%%%% the coefficient of forward waves in waveguend space
q = q_back;
p = p_back;





typeOfCylinder = 'Gyrotropic';
% typeOfCylinder = 'Isotropic';
% typeOfCylinder = 'PerfectConductivity';

rho = [0.12:0.1:50] *a_0/10;
% rho = [0.99:0.001:1.01] * 1*a_0;

%%%%%%%%%%%%%%%%%%%%%%% E-polarized beam %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  AE_0 = 1;
  Ez_inc    = AE_0.* (besselj(m, k_0.* rho * q_0)); 
  dEz_q_inc = AE_0.* k_0.* (q_0).*((besselj(m,   k_0.* rho* q_0) * (m))./ (k_0.* rho* q_0)  - besselj(m + 1,   k_0.* rho* q_0));

%%%%%%%%%%%%%%%%%%%%%%% H-polarized beam %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  AH_0 = 0;
  Hz_inc    = AH_0.* (besselj(m, k_0.* rho * q_0));
  dHz_q_inc = AH_0.* k_0.* (q_0).*((besselj(m,   k_0.* rho* q_0) * (m))./ (k_0.* rho* q_0)  - besselj(m + 1,   k_0.* rho* q_0));
  
 %%%%%%%%%%%%%%%%%%%%%%% other field components of the beam %%%%%%%%%%%%%%%
  A0 = 1./ (k_0 * (1 - p_0.^2));
  Ephi_inc  = A0.* ((-p_0.* (m./rho)).* Ez_inc + 1i * dHz_q_inc);
  Hphi_inc  = A0.* (-1i * dEz_q_inc - (p_0.* (m./rho)).* Hz_inc);
  Hrho_inc  = A0.* ((m./rho).* Ez_inc - 1i * p_0.* dHz_q_inc);
  Erho_inc  = A0.* (-1i * p_0.* dEz_q_inc - (m./rho).* Hz_inc);
  
%   plot_all_field_components(2, rho, a_0, [], Ez_inc, Ephi_inc, Erho_inc, Hz_inc, Hphi_inc, Hrho_inc)

%%%%%%%%%%%%%%%%%%%%%%% Calculation of coefficients of singular waves in
%%%%%%%%%%%%%%%%%%%%%%% cylinder and singular waves in free space generated
%%%%%%%%%%%%%%%%%%%%%%% from them

%%%%%% singular wave in waveguides
 H2m0  = besselh(m, 2, k_0.* a_0 * q_0);
 dH2m0 = (H2m0 * m)./ (k_0.* a_0 * q_0)  - besselh(m + 1, 2, k_0.* a_0 * q_0);
 H1m0  = besselh(m, 1, k_0.* a_0 * q_0);
 dH1m0 = (H1m0 * m)./ (k_0.* a_0 * q_0)  - besselh(m + 1, 1, k_0.* a_0 * q_0);
 [b_p_1_delta_forw, b_p_2_delta_forw, psi_delta_1, psi_delta_2,...
    B_1_delta_1,B_2_delta_1,Cm2_delta_1,  Dm2_delta_1,...
    B_1_delta_2,B_2_delta_2,Cm2_delta_2,  Dm2_delta_2] = coefficientsOf_deltaFunctionWave(typeOfCylinder, q_0, p_0, k_0,k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, AE_0, AH_0);

z=0;
Ez_waveguidespace   = field_ForBesselBeamRepresentation_deltaFunction('Ez',   typeOfCylinder, rho, q_0, q_0, k_0,k_0, a_0, EE1, GG1, HH1, MU1, EE2, MU2, EE, MU, m, z, AE_0, AH_0);
Ephi_waveguidespace = field_ForBesselBeamRepresentation_deltaFunction('Ephi', typeOfCylinder, rho, q_0, q_0, k_0,k_0, a_0, EE1, GG1, HH1, MU1, EE2, MU2, EE, MU, m, z, AE_0, AH_0);
Erho_waveguidespace = field_ForBesselBeamRepresentation_deltaFunction('Erho', typeOfCylinder, rho, q_0, q_0, k_0,k_0, a_0, EE1, GG1, HH1, MU1, EE2, MU2, EE, MU, m, z, AE_0, AH_0);
Hz_waveguidespace   = field_ForBesselBeamRepresentation_deltaFunction('Hz',   typeOfCylinder, rho, q_0, q_0, k_0,k_0, a_0, EE1, GG1, HH1, MU1, EE2, MU2, EE, MU, m, z, AE_0, AH_0);
Hphi_waveguidespace = field_ForBesselBeamRepresentation_deltaFunction('Hphi', typeOfCylinder, rho, q_0, q_0, k_0,k_0, a_0, EE1, GG1, HH1, MU1, EE2, MU2, EE, MU, m, z, AE_0, AH_0);
Hrho_waveguidespace = field_ForBesselBeamRepresentation_deltaFunction('Hrho', typeOfCylinder, rho, q_0, q_0, k_0,k_0, a_0, EE1, GG1, HH1, MU1, EE2, MU2, EE, MU, m, z, AE_0, AH_0);

plot_all_field_components(2, rho, a_0, 'o-', 1, Ez_waveguidespace, Ephi_waveguidespace, Erho_waveguidespace, Hz_waveguidespace, Hphi_waveguidespace, Hrho_waveguidespace)

Sz = real(-Ephi_waveguidespace.* conj(Hrho_waveguidespace) + Erho_waveguidespace.* conj(Hphi_waveguidespace));
figure(19); plot(rho, Sz, 'ro')

%%%%%% singular wave in free space
[a_p_field_Ewave_forw, a_p_field_Hwave_forw] = coefficientsOf_deltaFunctionWaveOfFreeSpace(typeOfCylinder, q_0, p_0,...
    k_0,k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m,...
    b_p_1_delta_forw, b_p_2_delta_forw, psi_delta_1, psi_delta_2,...
    Cm2_delta_1,  Dm2_delta_1,...
    Cm2_delta_2,  Dm2_delta_2);
%%%%%%%%%%%%%%%%%%%%%%% E-polarized beam %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Ez_freespace    = field_ofSingularWavesInFreeSpace_deltaFunction('Ez',    rho, q_0, p_0, k_0, m, z, a_p_field_Ewave_forw, a_p_field_Hwave_forw);
  Hz_freespace    = field_ofSingularWavesInFreeSpace_deltaFunction('Hz',    rho, q_0, p_0, k_0, m, z, a_p_field_Ewave_forw, a_p_field_Hwave_forw);
  Ephi_freespace  = field_ofSingularWavesInFreeSpace_deltaFunction('Ephi',  rho, q_0, p_0, k_0, m, z, a_p_field_Ewave_forw, a_p_field_Hwave_forw);
  Hphi_freespace  = field_ofSingularWavesInFreeSpace_deltaFunction('Hphi',  rho, q_0, p_0, k_0, m, z, a_p_field_Ewave_forw, a_p_field_Hwave_forw);
  Hrho_freespace  = field_ofSingularWavesInFreeSpace_deltaFunction('Hrho',  rho, q_0, p_0, k_0, m, z, a_p_field_Ewave_forw, a_p_field_Hwave_forw);
  Erho_freespace  = field_ofSingularWavesInFreeSpace_deltaFunction('Erho',  rho, q_0, p_0, k_0, m, z, a_p_field_Ewave_forw, a_p_field_Hwave_forw);

% plot_all_field_components(2, rho, a_0, '-', 2, Ez_freespace, Ephi_freespace, Erho_freespace, Hz_freespace, Hphi_freespace, Hrho_freespace)

    
%%%%%%%%%%%%%% propagation constants of eigenwaves %%%%%%%%%%%%%%%%%%%%%%%%
if(strcmp(typeOfCylinder, 'Isotropic'))
    p_n__of_descreteMode_of_isotropicCyl
elseif(strcmp(typeOfCylinder, 'Gyrotropic'))
    p_n__of_descreteMode_of_gyrotropicCyl
end

% p_n = [];
       
    q_n = sqrt(1-p_n.^2);
    q_n = q_n.* (2*(imag(q_n) <= 0)-1);
    

% %%%%%%%%%%%%%% calculation of the first approximation and the parameters for
% %%%%%%%%%%%%%% calculation of the field 
[b_p_field_1_forw, b_p_field_1_back, b_p_field_2_forw, b_p_field_2_back,...
          psi_forward_1, psi_backward_1, psi_forward_1_transp, psi_backward_1_transp,...
          psi_forward_2, psi_backward_2, psi_forward_2_transp, psi_backward_2_transp,...
          B_1_forward_1,B_2_forward_1,Cm2_forward_1,  Dm2_forward_1,...
          B_1_backward_1,B_2_backward_1,Cm2_backward_1, Dm2_backward_1,...
          B_1_forward_2,B_2_forward_2,Cm2_forward_2,  Dm2_forward_2,...
          B_1_backward_2,B_2_backward_2,Cm2_backward_2, Dm2_backward_2] =...
          coeffsOfFieldAndExcitationCoeffs_of_Continues_discreteRepr(typeOfCylinder, q, q_0, k_0,k_0, a_0, EE1, GG1, HH1, MU1, EE2, MU2, EE, MU, m, z, AE_0, AH_0);




%%%%%%%%%%%%%%%%%%%%%%%%%%% test - field representation of waveguide waves
%%%%%%%%%%%%%%%%%%%%%%%%%%% as waves of free space
z = -0*a_0;
[a_p_Ebeam_forw, a_p_Hbeam_forw, a_p_Ebeam_back, a_p_Hbeam_back] =...
               scatteringCoeffs_of_gyrotrop_forFieldRepresentation(typeOfCylinder, dq_simp,...
               q, p, q_n, p_n, q_0, p_0, k_0,k_0, a_0, EE1, GG1, HH1, MU1, EE2, MU2, EE, MU, m, z, AE_0, AH_0);
% %%%%%%%%%%% For free space part %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_p_Ebeam_freespace = a_p_Ebeam_forw;
a_p_Hbeam_freespace = a_p_Hbeam_forw;

figure(15);hold on;plot(q,a_p_Ebeam_forw, 'ro-');hold off
figure(17);hold on;plot(q,a_p_Hbeam_forw, 'ro-');hold off

a_p_Ebeam_freespace(isnan(a_p_Ebeam_freespace)) = 0;
a_p_Hbeam_freespace(isnan(a_p_Hbeam_freespace)) = 0;
    
% Ez_freespace   = zeros(size(rho));
% Ephi_freespace = zeros(size(rho));
% Erho_freespace = zeros(size(rho));
% Hz_freespace   = zeros(size(rho));
% Hphi_freespace = zeros(size(rho));
% Hrho_freespace = zeros(size(rho));
for ir = 1:size(rho,2)
    Ez_freespace(ir) = Ez_freespace(ir) +...
          (k_0/2/pi) * sum(dq_simp.* field__ForBeamDiffraction_InFreeSpace_discreteRepr('Ez', typeOfCylinder,rho(ir), z, q, p, k_0,  m,...
          a_p_Ebeam_freespace, a_p_Hbeam_freespace).* exp(-1i * k_0 * (p) * z));
      
    Hz_freespace(ir) = Hz_freespace(ir) +...
          (k_0/2/pi) * sum(dq_simp.* field__ForBeamDiffraction_InFreeSpace_discreteRepr('Hz', typeOfCylinder,rho(ir), z, q, p, k_0,  m,...
          a_p_Ebeam_freespace, a_p_Hbeam_freespace).* exp(-1i * k_0 * (p) * z));
      
    Ephi_freespace(ir) = Ephi_freespace(ir) +...
          (k_0/2/pi) * sum(dq_simp.* field__ForBeamDiffraction_InFreeSpace_discreteRepr('Ephi', typeOfCylinder,rho(ir), z, q, p, k_0,  m,...
          a_p_Ebeam_freespace, a_p_Hbeam_freespace).* exp(-1i * k_0 * (p) * z));
              
    Hphi_freespace(ir) = Hphi_freespace(ir) +...
          (k_0/2/pi) * sum(dq_simp.* field__ForBeamDiffraction_InFreeSpace_discreteRepr('Hphi', typeOfCylinder,rho(ir), z, q, p, k_0,  m,...
          a_p_Ebeam_freespace, a_p_Hbeam_freespace).* exp(-1i * k_0 * (p) * z));
      
    Erho_freespace(ir) = Erho_freespace(ir) +...
          (k_0/2/pi) * sum(dq_simp.* field__ForBeamDiffraction_InFreeSpace_discreteRepr('Erho', typeOfCylinder,rho(ir), z, q, p, k_0,  m,...
          a_p_Ebeam_freespace, a_p_Hbeam_freespace).* exp(-1i * k_0 * (p) * z));
      
    Hrho_freespace(ir) = Hrho_freespace(ir) +...
          (k_0/2/pi) * sum(dq_simp.* field__ForBeamDiffraction_InFreeSpace_discreteRepr('Hrho', typeOfCylinder,rho(ir), z, q, p, k_0,  m,...
          a_p_Ebeam_freespace, a_p_Hbeam_freespace).* exp(-1i * k_0 * (p) * z));
end



a_p_Ebeam_freespace = a_p_Ebeam_back;
a_p_Hbeam_freespace = a_p_Hbeam_back;

a_p_Ebeam_freespace(isnan(a_p_Ebeam_freespace)) = 0;
a_p_Hbeam_freespace(isnan(a_p_Hbeam_freespace)) = 0;
    
for ir = 1:size(rho,2)
    Ez_freespace(ir) = Ez_freespace(ir) +...
          (k_0/2/pi) * sum(dq_simp.* field__ForBeamDiffraction_InFreeSpace_discreteRepr('Ez', typeOfCylinder,rho(ir), z, q, -p, k_0,  m,...
          a_p_Ebeam_freespace, a_p_Hbeam_freespace).* exp(-1i * k_0 * (-p) * z));
      
    Hz_freespace(ir) = Hz_freespace(ir) +...
          (k_0/2/pi) * sum(dq_simp.* field__ForBeamDiffraction_InFreeSpace_discreteRepr('Hz', typeOfCylinder,rho(ir), z, q, -p, k_0,  m,...
          a_p_Ebeam_freespace, a_p_Hbeam_freespace).* exp(-1i * k_0 * (-p) * z));
      
    Ephi_freespace(ir) = Ephi_freespace(ir) +...
          (k_0/2/pi) * sum(dq_simp.* field__ForBeamDiffraction_InFreeSpace_discreteRepr('Ephi', typeOfCylinder,rho(ir), z, q, -p, k_0,  m,...
          a_p_Ebeam_freespace, a_p_Hbeam_freespace).* exp(-1i * k_0 * (-p) * z));
              
    Hphi_freespace(ir) = Hphi_freespace(ir) +...
          (k_0/2/pi) * sum(dq_simp.* field__ForBeamDiffraction_InFreeSpace_discreteRepr('Hphi', typeOfCylinder,rho(ir), z, q, -p, k_0,  m,...
          a_p_Ebeam_freespace, a_p_Hbeam_freespace).* exp(-1i * k_0 * (-p) * z));
      
    Erho_freespace(ir) = Erho_freespace(ir) +...
          (k_0/2/pi) * sum(dq_simp.* field__ForBeamDiffraction_InFreeSpace_discreteRepr('Erho', typeOfCylinder,rho(ir), z, q, -p, k_0,  m,...
          a_p_Ebeam_freespace, a_p_Hbeam_freespace).* exp(-1i * k_0 * (-p) * z));
      
    Hrho_freespace(ir) = Hrho_freespace(ir) +...
          (k_0/2/pi) * sum(dq_simp.* field__ForBeamDiffraction_InFreeSpace_discreteRepr('Hrho', typeOfCylinder,rho(ir), z, q, -p, k_0,  m,...
          a_p_Ebeam_freespace, a_p_Hbeam_freespace).* exp(-1i * k_0 * (-p) * z));
end

plot_all_field_components(2, rho, a_0, '-', 2, Ez_freespace, Ephi_freespace, Erho_freespace, Hz_freespace, Hphi_freespace, Hrho_freespace)


% figure(87)
% hold on
% plot(q, abs(a_p_Ebeam_forw))
% plot(q, abs(a_p_Hbeam_forw), 'r')
% hold off
% 
% figure(88)
% hold on
% plot(q, abs(a_p_Ebeam_back))
% plot(q, abs(a_p_Hbeam_back), 'r')
% hold off

Sz = real(-Ephi_freespace.* conj(Hrho_freespace) + Erho_freespace.* conj(Hphi_freespace));
figure(19); hold on;plot(rho, Sz,'b');hold off;


toc
  
  
