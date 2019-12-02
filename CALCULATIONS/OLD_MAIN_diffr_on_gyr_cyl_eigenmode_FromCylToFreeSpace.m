%%%%%% diffraction of bessel beams E and H polarized on edge of cylinder
%%%%%% edge is located at cross section z0=L 

clear all
tic
systemParameters
% EE1 = 1.002;
% HH1 = EE1;

L = 0;
z0 = 0;
%%%%%% paramters of beam
% m = 1;
p_0 = 1/sqrt(2);

p_0 = [
% 29.8516418826982;
31.3317352095852; %%% n=7, m=1, w_0/w_H = 0.025;
% 94.7689456855476 %%% n=35, m=1, w_0/w_H = 0.025;

% 94.5592839178522  %%% n=35, m=1, w_0/w_H = 0.022;
% 90.3451314635158  %%% n=35, m=1, w_0/w_H = 0.01;
% 87.0641528290118;
          ];
      
% p_0 = [
% %       3.08945;
% %       2.8356;
% %       2.7446;
%       2.2196;
% %       1.97
% %       1.010359864
%           ]; %%%a_0 = 4.330161354189217;k_0 = 0.541270169273652 *(2*pi)/a_0; EE1 = 10; m=1

q_0 = sqrt(1-p_0^2);
q_0 = q_0.* (2*(imag(q_0) <= 0)-1);

%%%%%% paramenter of cylinder for perfect conductivity core
% % a_0 = 2.4048./ (k_0 * q_0);%%% radius of inner core
% a_0 = 4.330161354189217;
% % k_0 = 0.541270169273652 * (2*pi) / a_0;
% % k_0 = 0.2 * (2*pi) / a_0;
% a_0 = 1.6;

q_0 = sqrt(1-p_0^2);
q_0 = q_0.* (2*(imag(q_0) <= 0)-1);

% q_0 = 0.5 + 1e-5;
% p_0 = sqrt(1 - q_0.^2);
% p_0 = real(p_0) - 1i * abs(imag(p_0));

typeOfCylinder = 'Gyrotropic';
% typeOfCylinder = 'Isotropic';
% typeOfCylinder = 'PerfectConductivity';

rho = [0.12:0.01:5] * a_0;

Exp_kp0z0 = exp(1i * k_0 * p_0 * (-L));

%%%%%%%%%%%%%%%%%%%%%%% E-polarized beam %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  AE_0 = 1.* Exp_kp0z0;
  Ez_inc    = AE_0.* (besselj(m, k_0.* rho * q_0)); 
  dEz_q_inc = AE_0.* k_0.* (q_0).* ((besselj(m,   k_0.* rho* q_0) * (m))./ (k_0.* rho* q_0)  - besselj(m + 1,   k_0.* rho* q_0));

%%%%%%%%%%%%%%%%%%%%%%% H-polarized beam %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  AH_0 = 0.* Exp_kp0z0;
  Hz_inc    = AH_0.* (besselj(m, k_0.* rho * q_0));
  dHz_q_inc = AH_0.* k_0.* (q_0).* ((besselj(m,   k_0.* rho* q_0) * (m))./ (k_0.* rho* q_0)  - besselj(m + 1,   k_0.* rho* q_0));
  
 %%%%%%%%%%%%%%%%%%%%%%% other field components of the beam %%%%%%%%%%%%%%%
  A0 = 1./ (k_0 * (1 - p_0.^2));
  Ephi_inc  = A0.* ((-p_0.* (m./rho)).* Ez_inc + 1i * dHz_q_inc);
  Hphi_inc  = A0.* (-1i * dEz_q_inc - (p_0.* (m./rho)).* Hz_inc);
  Hrho_inc  = A0.* ((m./rho).* Ez_inc - 1i * p_0.* dHz_q_inc);
  Erho_inc  = A0.* (-1i * p_0.* dEz_q_inc - (m./rho).* Hz_inc);
 
%   plot_all_field_components(3, rho, a_0, '-*', 1, Ez_inc, Ephi_inc, Erho_inc, Hz_inc, Hphi_inc, Hrho_inc)

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
    B_1_delta_2,B_2_delta_2,Cm2_delta_2,  Dm2_delta_2] = coefficientsOf_deltaFunctionWave(typeOfCylinder, q_0, p_0, k_0, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, AE_0, AH_0);

%%%%% ????? ????? z0 = 0, ?.?. ?????g???????? ?????? ????????? ?
%%%%% ???????????? ??????????? ?? ?????? ?????
z=0;
Ez_waveguidespace   = field_ForBesselBeamRepresentation_deltaFunction('Ez',   typeOfCylinder, rho, q_0, q_0, k_0,k_0, a_0, EE1, GG1, HH1, MU1, EE2, MU2, EE, MU, m, z, AE_0, AH_0);
Ephi_waveguidespace = field_ForBesselBeamRepresentation_deltaFunction('Ephi', typeOfCylinder, rho, q_0, q_0, k_0,k_0, a_0, EE1, GG1, HH1, MU1, EE2, MU2, EE, MU, m, z, AE_0, AH_0);
Erho_waveguidespace = field_ForBesselBeamRepresentation_deltaFunction('Erho', typeOfCylinder, rho, q_0, q_0, k_0,k_0, a_0, EE1, GG1, HH1, MU1, EE2, MU2, EE, MU, m, z, AE_0, AH_0);
Hz_waveguidespace   = field_ForBesselBeamRepresentation_deltaFunction('Hz',   typeOfCylinder, rho, q_0, q_0, k_0,k_0, a_0, EE1, GG1, HH1, MU1, EE2, MU2, EE, MU, m, z, AE_0, AH_0);
Hphi_waveguidespace = field_ForBesselBeamRepresentation_deltaFunction('Hphi', typeOfCylinder, rho, q_0, q_0, k_0,k_0, a_0, EE1, GG1, HH1, MU1, EE2, MU2, EE, MU, m, z, AE_0, AH_0);
Hrho_waveguidespace = field_ForBesselBeamRepresentation_deltaFunction('Hrho', typeOfCylinder, rho, q_0, q_0, k_0,k_0, a_0, EE1, GG1, HH1, MU1, EE2, MU2, EE, MU, m, z, AE_0, AH_0);

% plot_all_field_components(3, rho, a_0, 'o-', 1, Ez_waveguidespace, Ephi_waveguidespace, Erho_waveguidespace, Hz_waveguidespace, Hphi_waveguidespace, Hrho_waveguidespace)

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

% plot_all_field_components(3, rho, a_0, '-', 2, Ez_freespace, Ephi_freespace, Erho_freespace, Hz_freespace, Hphi_freespace, Hrho_freespace)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% calculation of scattering coefficients %%%%%%%%%%%%%%%%%%%
upper_Bound = 15;
dq = 0.01;
N = 50;

upper_Bound = 8000;
N = 800;

upper_Bound = 400;
N = 4;

% [dq_simp, q_back, p] = quadratureMethod_forIntegralEqs('trapeciech', N, q_0, 1e-8, upper_Bound);
[dq_simp, q_back, p] = quadratureMethod_forIntegralEqs('simpson', N, q_0, 1e-8, upper_Bound);

p_back = sqrt(1 - q_back.^2);
p_back = real(p_back) - 1i * abs(imag(p_back));
    
%%%%% the coefficient of forward waves in waveguend space
q = q_back;
p = p_back;

   
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
%       
% % % figure(1);hold on; plot(q, real(k_0*b_p_field_1_forw));hold off
% for in = 1:size(q_n,1)
%     b_n_forw(in,1) = a_p_field_descreteMode__ForBesselBeam(typeOfCylinder, q_n(in), q_0,  p_n(in), k_0, k_0, a_0, EE1, GG1, HH1, MU1, EE2, MU2, EE, MU, m, m, z, AE_0, AH_0);
%     b_n_back(in,1) = a_p_field_descreteMode__ForBesselBeam(typeOfCylinder, q_n(in), q_0, -p_n(in), k_0, k_0, a_0, EE1, GG1, HH1, MU1, EE2, MU2, EE, MU, m, m, z, AE_0, AH_0);
%     N_p_beta_des(in,1) = Norm_of_descreteMode(typeOfCylinder, q_n(in), p_n(in), k_0, a_0, EE1, GG1, HH1, MU1, EE2, MU2, EE, MU, m);
% end    
    
    
%%%%%%%%%%%%%%% scattering coefficients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% quadrature method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[b_p_field_1_back, b_p_field_2_back, b_n_back, a_p_Ebeam_forw,a_p_Hbeam_forw] =...           
          scatteringCoeffs_of_gyrotropCylinderEdge_quadMethod(typeOfCylinder, dq_simp,...
               q, p, q_n, p_n, q_0, p_0, k_0, k_0, a_0, EE1, GG1, HH1, MU1, EE2, MU2, EE, MU, m, 0, AE_0, AH_0);
           
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%% painting of field %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%% For space with waveguide %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 

%%%%%%%%%%%%%%%%%%% add delta function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ez_waveguidespace   = field_ForBesselBeamRepresentation_deltaFunction('Ez',   typeOfCylinder, rho, q_0, q_0, k_0,k_0, a_0, EE1, GG1, HH1, MU1, EE2, MU2, EE, MU, m, z, AE_0, AH_0);
Ephi_waveguidespace = field_ForBesselBeamRepresentation_deltaFunction('Ephi', typeOfCylinder, rho, q_0, q_0, k_0,k_0, a_0, EE1, GG1, HH1, MU1, EE2, MU2, EE, MU, m, z, AE_0, AH_0);
Erho_waveguidespace = field_ForBesselBeamRepresentation_deltaFunction('Erho', typeOfCylinder, rho, q_0, q_0, k_0,k_0, a_0, EE1, GG1, HH1, MU1, EE2, MU2, EE, MU, m, z, AE_0, AH_0);
Hz_waveguidespace   = field_ForBesselBeamRepresentation_deltaFunction('Hz',   typeOfCylinder, rho, q_0, q_0, k_0,k_0, a_0, EE1, GG1, HH1, MU1, EE2, MU2, EE, MU, m, z, AE_0, AH_0);
Hphi_waveguidespace = field_ForBesselBeamRepresentation_deltaFunction('Hphi', typeOfCylinder, rho, q_0, q_0, k_0,k_0, a_0, EE1, GG1, HH1, MU1, EE2, MU2, EE, MU, m, z, AE_0, AH_0);
Hrho_waveguidespace = field_ForBesselBeamRepresentation_deltaFunction('Hrho', typeOfCylinder, rho, q_0, q_0, k_0,k_0, a_0, EE1, GG1, HH1, MU1, EE2, MU2, EE, MU, m, z, AE_0, AH_0);

 for in=1:size(p_n,1)
        Ez_waveguidespace = Ez_waveguidespace +...
            field_of_discreteMode__ForBesselBeamRepresentation_OneDirection('Ez', typeOfCylinder, rho, q_n(in), -p_n(in), b_n_back(in),...
            q_0, k_0,k_0, a_0, EE1, GG1, HH1, MU1, EE2, MU2, EE, MU, m, z, AE_0, AH_0);
        Hz_waveguidespace = Hz_waveguidespace +...
            field_of_discreteMode__ForBesselBeamRepresentation_OneDirection('Hz', typeOfCylinder, rho, q_n(in), -p_n(in), b_n_back(in),...
            q_0, k_0,k_0, a_0, EE1, GG1, HH1, MU1, EE2, MU2, EE, MU, m, z, AE_0, AH_0);
        Ephi_waveguidespace = Ephi_waveguidespace +...
            field_of_discreteMode__ForBesselBeamRepresentation_OneDirection('Ephi', typeOfCylinder, rho, q_n(in), -p_n(in), b_n_back(in),...
            q_0, k_0,k_0, a_0, EE1, GG1, HH1, MU1, EE2, MU2, EE, MU, m, z, AE_0, AH_0);
        Erho_waveguidespace = Erho_waveguidespace +...
            field_of_discreteMode__ForBesselBeamRepresentation_OneDirection('Erho', typeOfCylinder, rho, q_n(in), -p_n(in), b_n_back(in),...
            q_0, k_0,k_0, a_0, EE1, GG1, HH1, MU1, EE2, MU2, EE, MU, m, z, AE_0, AH_0);  
        Hphi_waveguidespace = Hphi_waveguidespace +...
            field_of_discreteMode__ForBesselBeamRepresentation_OneDirection('Hphi', typeOfCylinder, rho, q_n(in), -p_n(in), b_n_back(in),...
            q_0, k_0,k_0, a_0, EE1, GG1, HH1, MU1, EE2, MU2, EE, MU, m, z, AE_0, AH_0);
        Hrho_waveguidespace = Hrho_waveguidespace +...
            field_of_discreteMode__ForBesselBeamRepresentation_OneDirection('Hrho', typeOfCylinder, rho, q_n(in), -p_n(in), b_n_back(in),...
            q_0, k_0,k_0, a_0, EE1, GG1, HH1, MU1, EE2, MU2, EE, MU, m, z, AE_0, AH_0);
 end

% b_p_field_1_forw = b_p_field_1_forw*0;
% b_p_field_2_forw = b_p_field_2_forw*0;
b_p_field_1_forw = zeros(1,size(b_p_field_1_forw,2));
b_p_field_2_forw = zeros(1,size(b_p_field_2_forw,2));
for ir = 1:size(rho,2)
    Ez_waveguidespace(ir) = Ez_waveguidespace(ir) +...
          (k_0/2/pi) * sum(dq_simp.* field__ForBesselBeamRepresentation_Continues_discreteRepr('Ez', typeOfCylinder,...
                  rho(ir), q, q_0, k_0,k_0, a_0, EE1, GG1, HH1, MU1, EE2, MU2, EE, MU, m, z, AE_0, AH_0,...
          b_p_field_1_forw, b_p_field_1_back, b_p_field_2_forw, b_p_field_2_back,...
          psi_forward_1, psi_backward_1, psi_forward_1_transp, psi_backward_1_transp,...
          psi_forward_2, psi_backward_2, psi_forward_2_transp, psi_backward_2_transp,...
          B_1_forward_1,B_2_forward_1,Cm2_forward_1,  Dm2_forward_1,...
          B_1_backward_1,B_2_backward_1,Cm2_backward_1, Dm2_backward_1,...
          B_1_forward_2,B_2_forward_2,Cm2_forward_2,  Dm2_forward_2,...
          B_1_backward_2,B_2_backward_2,Cm2_backward_2, Dm2_backward_2));
    Hz_waveguidespace(ir) = Hz_waveguidespace(ir) +...
          (k_0/2/pi) * sum(dq_simp.* field__ForBesselBeamRepresentation_Continues_discreteRepr('Hz', typeOfCylinder,...
                  rho(ir), q, q_0, k_0,k_0, a_0, EE1, GG1, HH1, MU1, EE2, MU2, EE, MU, m, z, AE_0, AH_0,...
          b_p_field_1_forw, b_p_field_1_back, b_p_field_2_forw, b_p_field_2_back,...
          psi_forward_1, psi_backward_1, psi_forward_1_transp, psi_backward_1_transp,...
          psi_forward_2, psi_backward_2, psi_forward_2_transp, psi_backward_2_transp,...
          B_1_forward_1,B_2_forward_1,Cm2_forward_1,  Dm2_forward_1,...
          B_1_backward_1,B_2_backward_1,Cm2_backward_1, Dm2_backward_1,...
          B_1_forward_2,B_2_forward_2,Cm2_forward_2,  Dm2_forward_2,...
          B_1_backward_2,B_2_backward_2,Cm2_backward_2, Dm2_backward_2));
    Ephi_waveguidespace(ir) = Ephi_waveguidespace(ir) +...
          (k_0/2/pi) * sum(dq_simp.* field__ForBesselBeamRepresentation_Continues_discreteRepr('Ephi', typeOfCylinder,...
                  rho(ir), q, q_0, k_0,k_0, a_0, EE1, GG1, HH1, MU1, EE2, MU2, EE, MU, m, z, AE_0, AH_0,...
          b_p_field_1_forw, b_p_field_1_back, b_p_field_2_forw, b_p_field_2_back,...
          psi_forward_1, psi_backward_1, psi_forward_1_transp, psi_backward_1_transp,...
          psi_forward_2, psi_backward_2, psi_forward_2_transp, psi_backward_2_transp,...
          B_1_forward_1,B_2_forward_1,Cm2_forward_1,  Dm2_forward_1,...
          B_1_backward_1,B_2_backward_1,Cm2_backward_1, Dm2_backward_1,...
          B_1_forward_2,B_2_forward_2,Cm2_forward_2,  Dm2_forward_2,...
          B_1_backward_2,B_2_backward_2,Cm2_backward_2, Dm2_backward_2));
              
    Hphi_waveguidespace(ir) = Hphi_waveguidespace(ir) +...
          (k_0/2/pi) * sum(dq_simp.* field__ForBesselBeamRepresentation_Continues_discreteRepr('Hphi', typeOfCylinder,...
                  rho(ir), q, q_0, k_0,k_0, a_0, EE1, GG1, HH1, MU1, EE2, MU2, EE, MU, m, z, AE_0, AH_0,...
          b_p_field_1_forw, b_p_field_1_back, b_p_field_2_forw, b_p_field_2_back,...
          psi_forward_1, psi_backward_1, psi_forward_1_transp, psi_backward_1_transp,...
          psi_forward_2, psi_backward_2, psi_forward_2_transp, psi_backward_2_transp,...
          B_1_forward_1,B_2_forward_1,Cm2_forward_1,  Dm2_forward_1,...
          B_1_backward_1,B_2_backward_1,Cm2_backward_1, Dm2_backward_1,...
          B_1_forward_2,B_2_forward_2,Cm2_forward_2,  Dm2_forward_2,...
          B_1_backward_2,B_2_backward_2,Cm2_backward_2, Dm2_backward_2));
    Erho_waveguidespace(ir) = Erho_waveguidespace(ir) +...
          (k_0/2/pi) * sum(dq_simp.* field__ForBesselBeamRepresentation_Continues_discreteRepr('Erho', typeOfCylinder,...
                  rho(ir), q, q_0, k_0,k_0, a_0, EE1, GG1, HH1, MU1, EE2, MU2, EE, MU, m, z, AE_0, AH_0,...
          b_p_field_1_forw, b_p_field_1_back, b_p_field_2_forw, b_p_field_2_back,...
          psi_forward_1, psi_backward_1, psi_forward_1_transp, psi_backward_1_transp,...
          psi_forward_2, psi_backward_2, psi_forward_2_transp, psi_backward_2_transp,...
          B_1_forward_1,B_2_forward_1,Cm2_forward_1,  Dm2_forward_1,...
          B_1_backward_1,B_2_backward_1,Cm2_backward_1, Dm2_backward_1,...
          B_1_forward_2,B_2_forward_2,Cm2_forward_2,  Dm2_forward_2,...
          B_1_backward_2,B_2_backward_2,Cm2_backward_2, Dm2_backward_2));
    Hrho_waveguidespace(ir) = Hrho_waveguidespace(ir)  +...
          (k_0/2/pi) * sum(dq_simp.* field__ForBesselBeamRepresentation_Continues_discreteRepr('Hrho', typeOfCylinder,...
                  rho(ir), q, q_0, k_0,k_0, a_0, EE1, GG1, HH1, MU1, EE2, MU2, EE, MU, m, z, AE_0, AH_0,...
          b_p_field_1_forw, b_p_field_1_back, b_p_field_2_forw, b_p_field_2_back,...
          psi_forward_1, psi_backward_1, psi_forward_1_transp, psi_backward_1_transp,...
          psi_forward_2, psi_backward_2, psi_forward_2_transp, psi_backward_2_transp,...
          B_1_forward_1,B_2_forward_1,Cm2_forward_1,  Dm2_forward_1,...
          B_1_backward_1,B_2_backward_1,Cm2_backward_1, Dm2_backward_1,...
          B_1_forward_2,B_2_forward_2,Cm2_forward_2,  Dm2_forward_2,...
          B_1_backward_2,B_2_backward_2,Cm2_backward_2, Dm2_backward_2));
end
plot_all_field_components(3, rho, a_0, 'o-', 1, Ez_waveguidespace, Ephi_waveguidespace, Erho_waveguidespace, Hz_waveguidespace, Hphi_waveguidespace, Hrho_waveguidespace)



% %%%%%%%%%%% For free space part %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_p_Ebeam_freespace = a_p_Ebeam_forw;
a_p_Hbeam_freespace = a_p_Hbeam_forw;

% a_p_Ebeam_freespace(isnan(a_p_Ebeam_freespace)) = 0;
% a_p_Hbeam_freespace(isnan(a_p_Hbeam_freespace)) = 0;
    
% Ez_freespace   = zeros(size(rho));
% Ephi_freespace = zeros(size(rho));
% Erho_freespace = zeros(size(rho));
% Hz_freespace   = zeros(size(rho));
% Hphi_freespace = zeros(size(rho));
% Hrho_freespace = zeros(size(rho));
for ir = 1:size(rho,2)
    Ez_freespace(ir) = Ez_freespace(ir) +...
          (k_0/2/pi) * sum(dq_simp.* field__ForBeamDiffraction_InFreeSpace_discreteRepr('Ez', typeOfCylinder,rho(ir), z, q, p, k_0,  m,...
          a_p_Ebeam_freespace, a_p_Hbeam_freespace));
      
    Hz_freespace(ir) = Hz_freespace(ir) +...
          (k_0/2/pi) * sum(dq_simp.* field__ForBeamDiffraction_InFreeSpace_discreteRepr('Hz', typeOfCylinder,rho(ir), z, q, p, k_0,  m,...
          a_p_Ebeam_freespace, a_p_Hbeam_freespace));
      
    Ephi_freespace(ir) = Ephi_freespace(ir) +...
          (k_0/2/pi) * sum(dq_simp.* field__ForBeamDiffraction_InFreeSpace_discreteRepr('Ephi', typeOfCylinder,rho(ir), z, q, p, k_0,  m,...
          a_p_Ebeam_freespace, a_p_Hbeam_freespace));
              
    Hphi_freespace(ir) = Hphi_freespace(ir) +...
          (k_0/2/pi) * sum(dq_simp.* field__ForBeamDiffraction_InFreeSpace_discreteRepr('Hphi', typeOfCylinder,rho(ir), z, q, p, k_0,  m,...
          a_p_Ebeam_freespace, a_p_Hbeam_freespace));
      
    Erho_freespace(ir) = Erho_freespace(ir) +...
          (k_0/2/pi) * sum(dq_simp.* field__ForBeamDiffraction_InFreeSpace_discreteRepr('Erho', typeOfCylinder,rho(ir), z, q, p, k_0,  m,...
          a_p_Ebeam_freespace, a_p_Hbeam_freespace));
      
    Hrho_freespace(ir) = Hrho_freespace(ir) +...
          (k_0/2/pi) * sum(dq_simp.* field__ForBeamDiffraction_InFreeSpace_discreteRepr('Hrho', typeOfCylinder,rho(ir), z, q, p, k_0,  m,...
          a_p_Ebeam_freespace, a_p_Hbeam_freespace));
end

plot_all_field_components(3, rho, a_0, '-', 2, Ez_freespace, Ephi_freespace, Erho_freespace, Hz_freespace, Hphi_freespace, Hrho_freespace)

% plot_all_field_components(3, rho, a_0, '-', 2, Ez_inc+Ez_freespace, Ephi_inc+Ephi_freespace,...
%     Erho_inc+Erho_freespace, Hz_inc+Hz_freespace, Hphi_inc+Hphi_freespace, Hrho_inc+Hrho_freespace)
% 

% figure(44);hold on;plot(q, real(b_p_field_1_back), q, imag(b_p_field_1_back), 'r');hold off;
% figure(45);hold on;plot(q, real(b_p_field_2_back), q, imag(b_p_field_2_back), 'r'); hold off;


figure(77)
hold on
plot(q, abs(b_p_field_1_back),'b')
plot(q, abs(b_p_field_2_back), 'r')
hold off

figure(78)
hold on
plot(q, abs(a_p_Ebeam_freespace),'b')
plot(q, abs(a_p_Hbeam_freespace), 'r')
hold off

toc




% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%% calculation of field in section %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% tic
% componentOfFieldSource = 'Ez';
% componentOfFieldSource = 'Ephi';
% % componentOfFieldSource = 'Erho';
% % componentOfFieldSource = 'Ex';
% % componentOfFieldSource = 'Ey';
% % componentOfFieldSource = 'Hz';
% componentOfFieldSource = 'Hphi';
% % componentOfFieldSource = 'Hrho';
% % componentOfFieldSource = 'Hy';
% 
% %%% находим поле в данной плоскости xOz с координатой y=const для z>0
%     y  = 1e-20;
% %%%% задаём границы в которых будем строить поле
%     bound_x = 5/2 * a_0;
%     bound_z = 15/2 * a_0;    
% 
%     step_x = 0.01 * bound_x;    
%     step_z = 0.01 * bound_z;
%     z_layers = [- bound_z:step_z: bound_z];
%     x_layers = [- bound_x :step_x: bound_x];
%     [x, z] = meshgrid(x_layers,z_layers);
%     
%     rho    = sqrt(x.^2 + y.^2);
%     phi1   = atan(abs(y)./ abs(x));
%     phi    = (phi1).* (x>0).* (y>0); 
%     phi    = phi + (pi - phi1).* (x<=0).* (y>0); 
%     phi    = phi + (pi + phi1).* (x<=0).* (y<=0); 
%     phi    = phi + (2*pi - phi1).* (x>0).* (y<=0);
%     
%     tic
%     compHE = zeros(size(z,1),size(z,2));
%     for iz = 1:size(z,1)
%             if(z(iz,1) <= 0) continue;end
%         for ix = 1:size(x,2)
%             compHE(iz,ix) = compHE(iz,ix)+...
%                 (k_0/2/pi) * sum(dq_simp.* field__ForBeamDiffraction_InFreeSpace_discreteRepr(componentOfFieldSource, typeOfCylinder,...
%                 rho(iz, ix), z(iz, ix), q, p, k_0,  m,...
%                 a_p_Ebeam_freespace, a_p_Hbeam_freespace)).* exp(- 1i * m * phi(iz,ix));
%         end
%     end    
%     
% %%% находим поле в данной плоскости xOz с координатой y=const для z>0
%     for iz = 1:size(z,1)
%         if(z(iz,1) > 0) continue;end
%         for ix = 1:size(x,2)
%             compHE(iz,ix) = compHE(iz,ix)+...
%                 (k_0/2/pi) * sum(dq_simp.* field_AdapterTo_ContinuesFieldFunction(componentOfFieldSource, typeOfCylinder,...
%                 rho(iz,ix),  z(iz, ix), q, -p, k_0, a_0, EE1, GG1, HH1, MU1, EE2, MU2, EE, MU, m,...
%                 b_p_field_1_back, psi_backward_1, B_1_backward_1, B_2_backward_1, Cm2_backward_1,  Dm2_backward_1)).* exp(- 1i * m * phi(iz,ix));
%             
%             compHE(iz,ix) = compHE(iz,ix)+...
%                 (k_0/2/pi) * sum(dq_simp.* field_AdapterTo_ContinuesFieldFunction(componentOfFieldSource, typeOfCylinder,...
%                 rho(iz,ix),  z(iz, ix), q, -p, k_0, a_0, EE1, GG1, HH1, MU1, EE2, MU2, EE, MU, m,...
%                 b_p_field_2_back, psi_backward_2, B_1_backward_2, B_2_backward_2, Cm2_backward_2,  Dm2_backward_2)).* exp(- 1i * m * phi(iz,ix));
%         end
%     end
%     %%% adding of discrete modes
%     for in=1:size(p_n,1)
%         compHE(z<=0) = compHE(z<=0) +...
%             field_of_discreteMode__ForBesselBeamRepresentation_OneDirection(componentOfFieldSource, typeOfCylinder, rho(z<=0), q_n(in), -p_n(in), b_n_back(in),...
%             q_0, k_0,k_0, a_0, EE1, GG1, HH1, MU1, EE2, MU2, EE, MU, m, z(z<=0), AE_0, AH_0).* exp(- 1i * m * phi(z<=0));
%     end
% %     %%%% adding of incident wave
%     compHE(z<=0) = compHE(z<=0) +...
%         k_0/2/pi * field_of_discreteMode__ForBesselBeamRepresentation_OneDirection(componentOfFieldSource, typeOfCylinder, rho(z<=0), q_0, p_0, 1,...
%         q_0, k_0,k_0, a_0, EE1, GG1, HH1, MU1, EE2, MU2, EE, MU, m, z(z<=0), AE_0, AH_0).* exp(- 1i * m * phi(z<=0)); %%% here k_0/2/pi - amplitude of icident wave (see field_ForBesselBeamRepresentation_deltaFunction)
% % compHE(z<=0) = compHE(z<=0) +...
% % field_ForBesselBeamRepresentation_deltaFunction(componentOfFieldSource,   typeOfCylinder,...
% % rho(z<=0), q_0, q_0, k_0,k_0, a_0, EE1, GG1, HH1, MU1, EE2, MU2, EE, MU, m, z(z<=0), AE_0, AH_0).* exp(- 1i * m * phi(z<=0));
%     
%     pcolorReImModule(z / a_0, 'z / a_0', x / a_0, 'x / a_0', compHE, componentOfFieldSource)

    
    

toc






  
  
