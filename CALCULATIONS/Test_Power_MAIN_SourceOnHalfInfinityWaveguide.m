clearvars
clc

tic
systemParameters
toOhms = 898755178736.818;

z0 = 0;
q_0 = sqrt(1-p_0^2);
q_0 = q_0.* (2*(imag(q_0) <= 0)-1);

typeOfCylinder = 'Gyrotropic';
% typeOfCylinder = 'Isotropic';

rho = [0.12:0.01:5] * 1*a_0;


%%%%%%%%%%%%%% propagation constants of eigenwaves %%%%%%%%%%%%%%%%%%%%%%%%
if(strcmp(typeOfCylinder, 'Isotropic'))
    p_n__of_descreteMode_of_isotropicCyl
elseif(strcmp(typeOfCylinder, 'Gyrotropic'))
    p_n__of_descreteMode_of_gyrotropicCyl
end
q_n = sqrt(1-p_n.^2);
q_n = q_n.* (2*(imag(q_n) <= 0)-1);

upper_Bound = 32000;
% upper_Bound = 16000;
% upper_Bound = 8000;
N = 1600;

% upper_Bound = 400;
% N = 40;

% [dq_simp, q_cs, p_cs] = quadratureMethod_forIntegralEqs('simpson', N, 1i, 1e-8, upper_Bound); % the method used before
%%[dq_simp, q_cs, p_cs] = quadratureMethod_forIntegralEqs('difficult_simpson', N, q_0, 1e-7, upper_Bound);


%%%%%% NEW INTEGRATION METHOD %%%%%%%%%%%%%%%%%%%%%%%%%%
N = 80;        %%% must be even
N_upper = 400; %%% must be even
upper_Bound = 32000;
q_n0 = [2; 20; 800; 16000];
[dq_simp, q_cs, p_cs] = quadratureMethod_forIntegralEqs_multipleIntervals('difficult_simpsonMarkovGildendurgExpWithCollisions', N, N_upper, q_n0, 1e-8, upper_Bound);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% the calculation of excitation coefficients of the discrete
%%%%%%%%%%%%%%% modes of source on surface of cylinder %%%%%%%%%%%%%%%%%%%%
for in = 1:size(q_n,1)
    a_smn_forw(in,1) = a_smn_of_discreteMode_Decorator(typeOfCylinder, q_n(in),  p_n(in), waveguideParameters, sourceParameters);
    a_smn_back(in,1) = a_smn_of_discreteMode_Decorator(typeOfCylinder, q_n(in), -p_n(in), waveguideParameters, sourceParameters);
end

for iq = 1:size(q_cs,2)
    a_cs_alfa1_forw(iq) = a_sma_of_continuousWaves_alpha1_Decorator(typeOfCylinder, q_cs(iq),  p_cs(iq), waveguideParameters, sourceParameters);
    a_cs_alfa1_back(iq) = a_sma_of_continuousWaves_alpha1_Decorator(typeOfCylinder, q_cs(iq), -p_cs(iq), waveguideParameters, sourceParameters);
    a_cs_alfa2_forw(iq) = a_sma_of_continuousWaves_alpha2_Decorator(typeOfCylinder, q_cs(iq),  p_cs(iq), waveguideParameters, sourceParameters);
    a_cs_alfa2_back(iq) = a_sma_of_continuousWaves_alpha2_Decorator(typeOfCylinder, q_cs(iq), -p_cs(iq), waveguideParameters, sourceParameters);
end

coeffsOfField = coeffsOfField_of_Continues_discreteRepr(typeOfCylinder, q_cs, p_cs, waveguideParameters, sourceParameters);

%%%%%%%%%%%%% test diffraction
% a_smn_forw(p_n ~= 31.331735209585200) = a_smn_forw(p_n ~= 31.331735209585200).* 0;
% a_smn_forw(p_n == 31.331735209585200) = 1;
% a_smn_forw = a_smn_forw * 0;

% a_cs_alfa1_forw = 0*a_cs_alfa1_forw;
% a_cs_alfa2_forw = 0*a_cs_alfa2_forw;
interalQ = q_cs > 0.999;
% a_cs_alfa1_forw(interalQ) = 0*a_cs_alfa1_forw(interalQ);
% a_cs_alfa2_forw(interalQ) = 0*a_cs_alfa2_forw(interalQ);

% a_cs_alfa1_forw(1) = 1;
a_cs_alfa1_forw(isinf(a_cs_alfa1_forw)) = isnan(a_cs_alfa1_forw(isinf(a_cs_alfa1_forw))) * 0;
a_cs_alfa2_forw(isinf(a_cs_alfa2_forw)) = isnan(a_cs_alfa2_forw(isinf(a_cs_alfa2_forw))) * 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_cs_alfa1_back(isinf(a_cs_alfa1_back)) = isnan(a_cs_alfa1_back(isinf(a_cs_alfa1_back))) * 0;
a_cs_alfa2_back(isinf(a_cs_alfa2_back)) = isnan(a_cs_alfa2_back(isinf(a_cs_alfa2_back))) * 0;
a_cs_alfa1_back_prop = a_cs_alfa1_back;
a_cs_alfa2_back_prop = a_cs_alfa2_back;
a_cs_alfa1_back_prop(interalQ) = 0*a_cs_alfa1_back_prop(interalQ);
a_cs_alfa2_back_prop(interalQ) = 0*a_cs_alfa2_back_prop(interalQ);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% forward excited waves %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z =  2*pi / k_0 * (0);
%L = sourceParameters.zCoordinate;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda_0 = 2 * pi / k_0; 
L_0 = 5*a_0; L_END = 2*lambda_0; NL = 7; %56 - ���� % ��������� ��� ������� 1120 ����� �� ���� ����� �����
% LL = L_0:(L_END - L_0)/NL:L_END;
LL = lambda_0;%[0.2 0.5 1.5]*lambda_0;
P_mod_back = zeros(size(LL,2),1);
P_cs_back  = zeros(size(LL,2),1);
P_cs_forw  = zeros(size(LL,2),1);
tic
for il = 1:size(LL,2) % new cycle on L created
    %%% only for upper-hybrid range! %%%%%%%%%%
%     clear p_n
%     p_n__of_descreteMode_of_gyrotropicCyl
%     q_n = sqrt(1-p_n.^2);
%     q_n = q_n.* (2*(imag(q_n) <= 0)-1);
    
    
    L = - LL(il);
    sourceParameters.zCoordinate = - LL(il);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
k_0 = sourceParameters.k_0;
% Ez_forwardExcitedWave = zeros(1,size(rho,2));
% Ephi_forwardExcitedWave = zeros(1,size(rho,2));
% Hphi_forwardExcitedWave = zeros(1,size(rho,2));
% %%% discrete mode
% for in = 1:size(q_n,1)
%     Ez_forwardExcitedWave = Ez_forwardExcitedWave +...
%         field_of_discreteMode_Decorator('Ez',   typeOfCylinder, rho, z, q_n(in), p_n(in), a_smn_forw(in), waveguideParameters, sourceParameters).*...
%         exp(1i * k_0 * p_n(in) * L);
% 
%     Ephi_forwardExcitedWave = Ephi_forwardExcitedWave +...
%         field_of_discreteMode_Decorator('Ephi',   typeOfCylinder, rho, z, q_n(in), p_n(in), a_smn_forw(in), waveguideParameters, sourceParameters).*...
%         exp(1i * k_0 * p_n(in) * L);
%     
%     Hphi_forwardExcitedWave = Hphi_forwardExcitedWave +...
%         field_of_discreteMode_Decorator('Hphi',   typeOfCylinder, rho, z, q_n(in), p_n(in), a_smn_forw(in), waveguideParameters, sourceParameters).*...
%         exp(1i * k_0 * p_n(in) * L);
% end


% %%% continues spectrume
% Expk0pcsL = exp(1i * k_0 * p_cs * L);
% a_cs_alfa1_forw_Exp = a_cs_alfa1_forw.* exp(1i * k_0 * p_cs * L);
% a_cs_alfa2_forw_Exp = a_cs_alfa2_forw.* exp(1i * k_0 * p_cs * L);
% for ir = 1:size(rho,2)
%     Ez_forwardExcitedWave(ir) = Ez_forwardExcitedWave(ir) +...
%           (k_0/2/pi) * sum(dq_simp.* field__ContinuesWaves_discreteRepr('Ez', typeOfCylinder, rho(ir), z, q_cs, p_cs,...
%           a_cs_alfa1_forw_Exp, a_cs_alfa1_forw_Exp*0, a_cs_alfa2_forw_Exp, a_cs_alfa2_forw_Exp*0,...
%           waveguideParameters, sourceParameters, coeffsOfField));
%       
%     Ephi_forwardExcitedWave(ir) = Ephi_forwardExcitedWave(ir) +...
%           (k_0/2/pi) * sum(dq_simp.* field__ContinuesWaves_discreteRepr('Ephi', typeOfCylinder, rho(ir), z, q_cs, p_cs,...
%           a_cs_alfa1_forw_Exp, a_cs_alfa1_forw_Exp*0, a_cs_alfa2_forw_Exp, a_cs_alfa2_forw_Exp*0,...
%           waveguideParameters, sourceParameters, coeffsOfField));
%       
%     Hphi_forwardExcitedWave(ir) = Hphi_forwardExcitedWave(ir) +...
%           (k_0/2/pi) * sum(dq_simp.* field__ContinuesWaves_discreteRepr('Hphi', typeOfCylinder, rho(ir), z, q_cs, p_cs,...
%           a_cs_alfa1_forw_Exp, a_cs_alfa1_forw_Exp*0, a_cs_alfa2_forw_Exp, a_cs_alfa2_forw_Exp*0,...
%           waveguideParameters, sourceParameters, coeffsOfField));
% end

% figure(1)
% hold on
% plot(rho / a_0, real(Ez_forwardExcitedWave),rho / a_0, imag(Ez_forwardExcitedWave),'r')
% title('E_{z}');
% hold off
% figure(2)
% hold on
% plot(rho / a_0, real(Ephi_forwardExcitedWave),rho / a_0, imag(Ephi_forwardExcitedWave),'r')
% title('E_{\phi}');
% hold off
% figure(3)
% hold on
% plot(rho / a_0, real(Hphi_forwardExcitedWave),rho / a_0, imag(Hphi_forwardExcitedWave),'r')
% title('H_{\phi}');
% hold off



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% calculation of diffarcted waves %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


continuousSperctrumCharacters.dq_simp = dq_simp;
continuousSperctrumCharacters.q = q_cs;
continuousSperctrumCharacters.p = p_cs;
continuousSperctrumCharacters.a_cs_alfa1_forw = a_cs_alfa1_forw;
continuousSperctrumCharacters.a_cs_alfa2_forw = a_cs_alfa2_forw;

discreteSperctrumCharacters.q_n = q_n;
discreteSperctrumCharacters.p_n = p_n;
discreteSperctrumCharacters.a_smn_forw = a_smn_forw;

[b_p_field_1_back, b_p_field_2_back, b_n_back, a_p_Ebeam_forw,a_p_Hbeam_forw] =...           
          scatteringCoeffs_of_EigenWavesOfHalfInfinityGyrotropCylEdge(typeOfCylinder,...
               continuousSperctrumCharacters, discreteSperctrumCharacters,waveguideParameters, sourceParameters);

           

 

%%%%%%%% CALCULATING OF PARTIAL POWERS %%%%%%%%%%% added by Oleg
%%%%%%%% Ostafiychuk 14.01.2020
%clc

b_p_field_1_back_prop = [b_p_field_1_back(~interalQ) 0*b_p_field_1_back(interalQ)];
b_p_field_2_back_prop = [b_p_field_2_back(~interalQ) 0*b_p_field_2_back(interalQ)];

a_p_Ebeam_forw_prop = [a_p_Ebeam_forw(~interalQ) 0*a_p_Ebeam_forw(interalQ)];
a_p_Hbeam_forw_prop = [a_p_Hbeam_forw(~interalQ) 0*a_p_Hbeam_forw(interalQ)];

%%% this is only for UH range !!! %%%%%%%%%%%%%%%%%%%%%%
% clear p_n
% p_n_set_without_complex_modes_for_UH
% q_n = sqrt(1-p_n.^2);
% q_n = q_n.* (2*(imag(q_n) <= 0)-1);
% for iin = 1:size(q_n,1);
%     a_smn_back_prop(iin) = a_smn_back(iin + 6);
%     b_n_back_prop(iin) = b_n_back(iin+6);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% in negative direction of z-axis%%%%%%%%%%
P_n_minus = zeros(size(q_n,1),1);
for in = 1:size(q_n,1)
    P_n_minus(in) = P_of_discreteMode(typeOfCylinder, q_n(in), - p_n(in), k_0, a_0, ...
                                                EE1, GG1, HH1, MU1, EE, MU, m, a_smn_back(in).*exp(1i * k_0 * p_n(in) * abs(L)) + b_n_back(in));
                                            
    
end

P_mod_back(il) = - sum(P_n_minus);

Phase_exp_cs = exp(1i * k_0 * p_cs * abs(L));
Phase_exp_cs(interalQ) = 0;
P_cs_minus_1 = ...
      sum(dq_simp.*P_of_continuousWaves_alpha1(typeOfCylinder, q_cs, p_0, - p_cs, k_0, a_0, ...
                                 EE1, GG1, HH1, MU1, EE, MU, c, m, z, j_f, j_z, d, a_cs_alfa1_back_prop.* Phase_exp_cs  + b_p_field_1_back_prop));
P_cs_minus_2 = ...
      sum(dq_simp.*P_of_continuousWaves_alpha2(typeOfCylinder, q_cs, p_0, - p_cs, k_0, a_0, ...
                                 EE1, GG1, HH1, MU1, EE, MU, c, m, z, j_f, j_z, d, a_cs_alfa2_back_prop.* Phase_exp_cs + b_p_field_2_back_prop));

P_cs_back(il) = - P_cs_minus_1 - P_cs_minus_2;

%%%%%%%% in half-space z>0%%%%%%%%%%

%%% the norms of eigenwaves of free space %%%
% N_p1 = @(q, p) (-1)^(m+1) * c * p./ (k_0.^2 * q);
% N_p2 = @(q, p) (-1)^(m+2) * c * p./ (k_0.^2 * q);
% P_cs_plus_1 = 1/4.*sum(dq_simp.*abs(a_p_Ebeam_forw_prop).^2.*N_p1(q_cs, p_cs));
% P_cs_plus_2 = 1/4.*sum(dq_simp.*abs(a_p_Hbeam_forw_prop).^2.*N_p2(q_cs, p_cs));

P_p1 = @(q, p) 1/4.*(-1)^(m+1) * c * p./ (k_0.^2 * q);
P_p2 = @(q, p) 1/4.*(-1)^(m+2) * c * p./ (k_0.^2 * q)*(-1);
P_cs_plus_1 = sum(dq_simp.*abs(a_p_Ebeam_forw_prop).^2.*P_p1(q_cs, p_cs));
P_cs_plus_2 = sum(dq_simp.*abs(a_p_Hbeam_forw_prop).^2.*P_p2(q_cs, p_cs));
                             
P_cs_forw(il) = P_cs_plus_1 + P_cs_plus_2;   

%P_ratio = P_mod_minus/(P_cs_minus + P_cs_plus)

end

toc

%% calculation of source resistance in free space
delta = 0;
pp = @(q) sqrt(1 - q.^2);
% N_p1 = @(q, p) (-1)^(m+1) * c * p./ (k_0.^2 * q);
% N_p2 = @(q, p) (-1)^(m+2) * c * p./ (k_0.^2 * q);

%%%%% for free space
P_Ewave_plus  =  quadgk(@(q) abs(a_sma1_freeSpace(q, pp(q), k_0, p_0, a_0, m, c, j_z, j_f, d)).^2.*...
    P_p1(q, pp(q)), delta, 1);
P_Ewave_minus  = quadgk(@(q) abs(a_sma1_freeSpace(q, - pp(q), k_0, p_0, a_0, m, c, j_z, j_f, d)).^2.*...
    P_p1(q,-pp(q)), delta, 1);

P_Hwave_plus  =  quadgk(@(q) abs(a_sma2_freeSpace(q, pp(q), k_0, p_0, a_0, m, c, j_z, j_f, d)).^2.*...
    P_p2(q, pp(q)), delta, 1);
P_Hwave_minus  =  quadgk(@(q) abs(a_sma2_freeSpace(q, - pp(q), k_0, p_0, a_0, m, c, j_z, j_f, d)).^2.*...
    P_p2(q,- pp(q)), delta, 1);

P_free_space_E = P_Ewave_plus - P_Ewave_minus;
P_free_space_H = P_Hwave_plus - P_Hwave_minus;
P_free_space = P_free_space_E + P_free_space_H
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(7)
% plot(LL/lambda_0, P_mod_back/P_free_space, 'r'); hold on
% plot(LL/lambda_0, 5*P_cs_back/P_free_space, 'b'); hold on
% plot(LL/lambda_0, 5*P_cs_forw/P_free_space, 'k'); hold on

% figure(10)
% bar(p_n, -P_n_minus/P_free_space);

P_Ratio = P_mod_back/(P_cs_back + P_cs_forw);

P_mod_back_on_Pfreespace = P_mod_back / P_free_space
P_cs_forw_on_Pfreespace  = P_cs_forw / P_free_space
P_cs_back_on_Pfreespace  = P_cs_back / P_free_space










%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%  test of norm of the modes with using the integrals %%%%%%%%%%
tic 

a_b = 100 * a_0;
% %%%%%% NEW INTEGRATION METHOD %%%%%%%%%%%%%%%%%%%%%%%%%%
% N = 40;        %%% must be even
% N_upper = 40; %%% must be even
% upper_Bound = 32000;
% q_n0 = [2; 20; 800; 16000];
% [dq_simp, q_cs, p_cs] = quadratureMethod_forIntegralEqs_multipleIntervals('difficult_simpsonMarkovGildendurgExpWithCollisions', N, N_upper, q_n0, 1e-8, upper_Bound);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q = q_cs;
p = p_cs;

P_n = zeros(size(p_n,1), size(LL,2));
P_n_incidence = P_n;
P_n_scatted = P_n;
P_beta1 = zeros(size(LL,2),1);
P_beta2 = P_beta1;
P_E = P_beta1;
P_H = P_beta1;
Sum_Modes = P_beta1;



        %%% first type of waves
        psi_forward_1   = psi1_q__gyrotropic(k_0, k_0, a_0, EE1, GG1, HH1, m,   p, q);
        psi_backward_1  = psi1_q__gyrotropic(k_0, k_0, a_0, EE1, GG1, HH1, m,  -p, q);
        psi_forward_1_transp   = psi1_q__gyrotropic(k_0, k_0, a_0, EE1, -GG1, HH1, -m,   p, q);
        psi_backward_1_transp  = psi1_q__gyrotropic(k_0, k_0, a_0, EE1, -GG1, HH1, -m,  -p, q);
        [B_1_forward_1,B_2_forward_1,Cm2_forward_1,  Dm2_forward_1]  = coefficientsOfContinuousSpectrum(typeOfCylinder,k_0, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m,   p, q,   psi_forward_1, 1);
        [B_1_backward_1,B_2_backward_1,Cm2_backward_1, Dm2_backward_1] = coefficientsOfContinuousSpectrum(typeOfCylinder,k_0, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, (-p), q, psi_backward_1, 1);
        %%% second type of waves
        psi_forward_2  = psi2_q__gyrotropic(k_0, k_0, a_0, EE1, GG1, HH1, m,   p, q);
        psi_backward_2  = psi2_q__gyrotropic(k_0, k_0, a_0, EE1, GG1, HH1, m,  -p, q);
        psi_forward_2_transp   = psi2_q__gyrotropic(k_0, k_0, a_0, EE1, -GG1, HH1, -m,   p, q);
        psi_backward_2_transp  = psi2_q__gyrotropic(k_0, k_0, a_0, EE1, -GG1, HH1, -m,  -p, q);
       [B_1_forward_2,B_2_forward_2,Cm2_forward_2,  Dm2_forward_2]  =   coefficientsOfContinuousSpectrum_forPsi2(typeOfCylinder,k_0, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m,   p,  q,  psi_forward_2, 1);
        [B_1_backward_2,B_2_backward_2,Cm2_backward_2, Dm2_backward_2] = coefficientsOfContinuousSpectrum_forPsi2(typeOfCylinder,k_0, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, (-p), q, psi_backward_2, 1);





%%%%%%%%% numerical calculation of radiation power in free space %%%%%%%%%%
% psi_forward_1 = psi_forward_1(q_cs<1);
% psi_forward_2 = psi_forward_2(q_cs<1);
% b_p_field_1_back = b_p_field_1_back(q_cs<1);
% b_p_field_2_back = b_p_field_2_back(q_cs<1);
% a_p_Ebeam_forw = a_p_Ebeam_forw(q_cs<1);
% a_p_Hbeam_forw = a_p_Hbeam_forw(q_cs<1);
% p = p;
% dq_simp = dq_simp;
% q = q_cs;

upper_Bound = 5 * a_0;
N = 1000;
[dr_simp_for_n, rho_for_n] = quadratureMethod_forPower('simpson', N, 1e-8, upper_Bound);

upper_Bound = 900000000 * a_b;
N = 1000;
            [dr_simp, rho] = quadratureMethod_forPower('simpson', N, 1e-8, upper_Bound);

% upper_Bound = 9000000 * a_b;
N = 1000;
  [dr_simp_beta, rho_beta] = quadratureMethod_forPower('simpson', N, 1e-8, upper_Bound);



for il = 1:size(LL,2) % new cycle on L created
    
    Phase_exp_cs = exp(1i * k_0 * p_cs * abs(L));
    Phase_exp_cs(interalQ) = 0;
    
%%%% powers of transversed eigenwaves
     S_component1 = zeros(1,size(rho_for_n,2));
     S_component2 = zeros(1,size(rho_for_n,2));
     S_component3 = zeros(1,size(rho_for_n,2));
     [B_1,B_2,Cm2, Dm2]  = coefficientsOfContinuousSpectrum(typeOfCylinder, k_0, k_0, a_0, EE1,GG1,HH1, 1, 1, 1, m, p_n, q_n, 0,  1);
     b_n_back_dop = a_smn_back.*exp(1i * k_0 * p_n * abs(L)) + b_n_back;
     for in = 1:size(q_n,1)
        S_component1 = componentOfPointingVector_Waveguide('S_z', typeOfCylinder, rho_for_n, 0, 0, q_n(in), -p_n(in), 0, b_n_back_dop(in),...
           k_0, a_0, EE1, GG1, HH1, 1, 1, 1, 1, 1, m, c,...
           B_1(in),B_2(in),Cm2(in),  Dm2(in), 0);
        P_n_f(in) = 2*pi * sum(dr_simp_for_n.* rho_for_n.* S_component1);
                
        S_component2 = componentOfPointingVector_Waveguide('S_z', typeOfCylinder, rho_for_n, 0, 0, q_n(in), -p_n(in), 0,  b_n_back(in),...
           k_0, a_0, EE1, GG1, HH1, 1, 1, 1, 1, 1, m, c,...
           B_1(in),B_2(in),Cm2(in),  Dm2(in), 0);
        P_n_f_scat(in) = 2*pi * sum(dr_simp_for_n.* rho_for_n.* S_component2);
        
        S_component3 = componentOfPointingVector_Waveguide('S_z', typeOfCylinder, rho_for_n, 0, 0, q_n(in), p_n(in), 0,  a_smn_forw(in),...
           k_0, a_0, EE1, GG1, HH1, 1, 1, 1, 1, 1, m, c,...
           B_1(in),B_2(in),Cm2(in),  Dm2(in), 0);
        P_n_f_inc(in) = 2*pi * sum(dr_simp_for_n.* rho_for_n.* S_component3);
     end
     P_n(:,il) = P_n_f / P_free_space;
     P_n_incidence(:,il) = P_n_f_inc/ P_free_space;
     P_n_scatted(:,il)   = P_n_f_scat/ P_free_space;
     Sum_Modes(il) = sum(P_n_f) / P_free_space
     Sum_Modes_inc(il) = sum(P_n_f_inc) / P_free_space
     Sum_Modes_scat(il) = sum(P_n_f_scat) / P_free_space
     

     Sum_Modes_inc(il) - abs(Sum_Modes_scat(il))
     figure(5)
     plot(p_n(:,1), P_n_f_inc(1,:), p_n(:,1), -P_n_f_scat(1,:), 'r')
     


    %%% first type of waves
    b_p_field_1_back = a_cs_alfa1_back_prop.* Phase_exp_cs  + b_p_field_1_back_prop;
    S_component2 = zeros(1,size(rho_beta,2));
    [B_1,B_2,Cm2,  Dm2]  = coefficientsOfContinuousSpectrum(typeOfCylinder,k_0, k_0, a_0, EE1,GG1,HH1, 1, 1, 1, m, -p, q, psi_forward_1,  1);
    for ir = 1:size(rho_beta,2)
        S_component2(ir) = componentOfPointingVector_Waveguide('S_z', typeOfCylinder, rho_beta(ir), 0, 0, q, -p, psi_forward_1, b_p_field_1_back,...
            k_0, a_0, EE1,GG1,HH1, 1, 1, 1, 1, 1, m, c,...
            B_1,B_2,Cm2,  Dm2, dq_simp);
    end
    P_beta1(il)   = 2*pi * (sum(dr_simp_beta.* rho_beta.* S_component2)) / P_free_space
    
%     f = @(x) componentOfPointingVector_Waveguide('S_z', typeOfCylinder, x, 0, 0, q, -p, psi_forward_1, b_p_field_1_back,...
%             k_0, a_0, EE1,GG1,HH1, 1, 1, 1, 1, 1, m, c,...
%             B_1,B_2,Cm2, Dm2, dq_simp);
%     quadgk(f, 1e-8, 1, 'abstol',1e-7) / P_free_space
    
    %%% second type of waves
    b_p_field_2_back = a_cs_alfa2_back_prop.* Phase_exp_cs + b_p_field_2_back_prop;
    S_component3 = zeros(1,size(rho_beta,2),1);
    [B_1,B_2,Cm2,  Dm2]  = coefficientsOfContinuousSpectrum_forPsi2(typeOfCylinder,k_0, k_0, a_0, EE1,GG1,HH1, 1, 1, 1, m, -p, q, psi_forward_2,  1);
    for ir = 1:size(rho_beta,2)
        S_component3(ir) = componentOfPointingVector_Waveguide('S_z', typeOfCylinder, rho_beta(ir), 0, 0, q, -p, psi_forward_2, b_p_field_2_back,...
            k_0, a_0, EE1,GG1,HH1, 1, 1, 1, 1, 1, m, c,...
            B_1,B_2,Cm2,  Dm2, dq_simp);
    end
    P_beta2(il) = 2*pi * (sum(dr_simp_beta.* rho_beta.* S_component3)) / P_free_space
    
    %%%% powers of waves which reflected to free space
    %%% E type of waves
    S_component4 = zeros(1,size(rho_beta,2),1);
    for ir = 1:size(rho,2)
        S_component4(ir) = componentOfPointingVector_freeSpace('S_z', typeOfCylinder, rho(ir), 0, 0, q, p, a_p_Ebeam_forw_prop, 0,...
            k_0, m, c, dq_simp);
    end
    P_E(il)   = 2*pi * (sum(dr_simp.* rho.* S_component4)) / P_free_space
    
    %%% H - type of waves
    S_component5 = zeros(1,size(rho_beta,2),1);
    for ir = 1:size(rho,2)
        S_component5(ir) = componentOfPointingVector_freeSpace('S_z', typeOfCylinder, rho(ir), 0, 0, q, p, 0, a_p_Hbeam_forw_prop,...
            k_0, m, c, dq_simp);
    end
    P_H(il) = 2*pi * (sum(dr_simp.* rho.* S_component5)) / P_free_space
    
    totalPower(il) = Sum_Modes(il) + P_beta1(il) + P_beta2(il) + P_E(il) + P_H(il)
    
    abs(P_beta1(il)) + abs(P_beta2(il)) + abs(P_E(il)) + abs(P_H(il))

end
toc




%%%%%%%%%%%% end test the norms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%













% %%%%% for waveguide
% P_alp1_forw  =  quadgk(@(q) P_of_continuousWaves_alpha1(typeOfCylinder,...
%     q, p_0, pp(q), k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, 0, j_f, j_z, d), delta, 1);
% P_alp1_back  = -quadgk(@(q) P_of_continuousWaves_alpha1(typeOfCylinder,...
%     q, p_0,-pp(q), k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, 0, j_f, j_z, d), delta, 1);
% 
% P_alp2_forw  = -quadgk(@(q) P_of_continuousWaves_alpha2(typeOfCylinder,...
%     q, p_0, pp(q), k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, 0, j_f, j_z, d), delta, 1);
% P_alp2_back  = quadgk(@(q) P_of_continuousWaves_alpha2(typeOfCylinder,...
%     q, p_0,-pp(q), k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, 0, j_f, j_z, d), delta, 1);


% R_E_freespace = abs((P_Ewave_plus + P_Ewave_minus)) * 2 / I_0^2 * toOhms
% R_H_freespace = abs((P_Hwave_plus + P_Hwave_minus)) * 2 / I_0^2 * toOhms
% R_alp1        = abs((P_alp1_forw + P_alp1_back))  * 2 / I_0^2   * toOhms
% R_alp2        = abs((P_alp2_forw + P_alp2_back))  * 2 / I_0^2   * toOhms

% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
% 
% Ez_waveguidespace = Ez_forwardExcitedWave;
% Ephi_waveguidespace = Ephi_forwardExcitedWave;
% Hphi_waveguidespace = Hphi_forwardExcitedWave;
% 
% %%% add field field of diffracted waves - discrete spectrum waves
% for in = 1:size(q_n,1)
%     Ez_waveguidespace = Ez_waveguidespace +...
%         field_of_discreteMode_Decorator('Ez',   typeOfCylinder, rho, z, q_n(in), -p_n(in), b_n_back(in), waveguideParameters, sourceParameters);
% 
%     Ephi_waveguidespace = Ephi_waveguidespace +...
%         field_of_discreteMode_Decorator('Ephi',   typeOfCylinder, rho, z, q_n(in), -p_n(in), b_n_back(in), waveguideParameters, sourceParameters);
%     
%     Hphi_waveguidespace = Hphi_waveguidespace +...
%         field_of_discreteMode_Decorator('Hphi',   typeOfCylinder, rho, z, q_n(in), -p_n(in), b_n_back(in), waveguideParameters, sourceParameters);
% end   
% 
% %%% add field field of diffracted waves - continuous spectrum waves
% for ir = 1:size(rho,2)
%     Ez_waveguidespace(ir) = Ez_waveguidespace(ir) +...
%           (k_0/2/pi) * sum(dq_simp.* field__ContinuesWaves_discreteRepr('Ez', typeOfCylinder, rho(ir), z, q_cs, p_cs,...
%           b_p_field_1_back*0, b_p_field_1_back, b_p_field_2_back*0, b_p_field_2_back,...
%           waveguideParameters, sourceParameters, coeffsOfField));
%       
%     Ephi_waveguidespace(ir) = Ephi_waveguidespace(ir) +...
%           (k_0/2/pi) * sum(dq_simp.* field__ContinuesWaves_discreteRepr('Ephi', typeOfCylinder, rho(ir), z, q_cs, p_cs,...
%           b_p_field_1_back*0, b_p_field_1_back, b_p_field_2_back*0, b_p_field_2_back,...
%           waveguideParameters, sourceParameters, coeffsOfField));
%       
%     Hphi_waveguidespace(ir) = Hphi_waveguidespace(ir) +...
%           (k_0/2/pi) * sum(dq_simp.* field__ContinuesWaves_discreteRepr('Hphi', typeOfCylinder, rho(ir), z, q_cs, p_cs,...
%           b_p_field_1_back*0, b_p_field_1_back, b_p_field_2_back*0, b_p_field_2_back,...
%           waveguideParameters, sourceParameters, coeffsOfField));
% end
% 
% figure(1)
% hold on
% plot(rho / a_0, real(Ez_waveguidespace), 'b-o', rho / a_0, imag(Ez_waveguidespace),'r-o')
% title('E_{z}');
% hold off
% figure(2)
% hold on
% plot(rho / a_0, real(Ephi_waveguidespace), 'b-o', rho / a_0, imag(Ephi_waveguidespace),'r-o')
% title('E_{\phi}');
% hold off
% figure(3)
% hold on
% plot(rho / a_0, real(Hphi_waveguidespace), 'b-o', rho / a_0, imag(Hphi_waveguidespace),'r-o')
% title('H_{\phi}');
% hold off
% 
% % %%%%%%%%%%% For free space part %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a_p_Ebeam_freespace = a_p_Ebeam_forw;
% a_p_Hbeam_freespace = a_p_Hbeam_forw;
%     
% Ez_freespace   = zeros(size(rho));
% Ephi_freespace = zeros(size(rho));
% % Erho_freespace = zeros(size(rho));
% % Hz_freespace   = zeros(size(rho));
% Hphi_freespace = zeros(size(rho));
% % Hrho_freespace = zeros(size(rho));
% for ir = 1:size(rho,2)
%     Ez_freespace(ir) = Ez_freespace(ir) +...
%           (k_0/2/pi) * sum(dq_simp.* field__ForBeamDiffraction_InFreeSpace_discreteRepr('Ez', typeOfCylinder,rho(ir), z, q_cs, p_cs, k_0,  m,...
%           a_p_Ebeam_freespace, a_p_Hbeam_freespace));
%       
%     Ephi_freespace(ir) = Ephi_freespace(ir) +...
%           (k_0/2/pi) * sum(dq_simp.* field__ForBeamDiffraction_InFreeSpace_discreteRepr('Ephi', typeOfCylinder,rho(ir), z, q_cs, p_cs, k_0,  m,...
%           a_p_Ebeam_freespace, a_p_Hbeam_freespace));
%               
%     Hphi_freespace(ir) = Hphi_freespace(ir) +...
%           (k_0/2/pi) * sum(dq_simp.* field__ForBeamDiffraction_InFreeSpace_discreteRepr('Hphi', typeOfCylinder,rho(ir), z, q_cs, p_cs, k_0,  m,...
%           a_p_Ebeam_freespace, a_p_Hbeam_freespace));
% end
% 
% figure(1)
% hold on
% plot(rho / a_0, real(Ez_freespace), 'b', rho / a_0, imag(Ez_freespace),'r')
% title('E_{z}');
% hold off
% figure(2)
% hold on
% plot(rho / a_0, real(Ephi_freespace), 'b',rho / a_0, imag(Ephi_freespace),'r')
% title('E_{\phi}');
% hold off
% figure(3)
% hold on
% plot(rho / a_0, real(Hphi_freespace), 'b',rho / a_0, imag(Hphi_freespace),'r')
% title('H_{\phi}');
% hold off

toc
  
  
