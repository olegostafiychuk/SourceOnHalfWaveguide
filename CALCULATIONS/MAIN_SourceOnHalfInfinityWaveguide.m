clearvars
% clc

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
%     p_n = p_n(1:40);
    p_n(7:end) = -p_n(7:end);
end
[q1, q2, n1, n2, alp1, alp2, bet1, bet2] = term_of_gyrotropic_waveguide(EE1, GG1, HH1, p_n);
% figure(111);
% plot(real(p_n), real(q1), 'b*', real(p_n), real(q2),'r*')

% p_n_losses;%этот скрипт подт€гивает значени€ продольных волновых чисел и столбец столкновений из файла disperson_losses_m1.mat
% p_n_los = p_n_real + 1i*p_n_imag;
% p_n = p_n_los(1:240,71);
% Nu_e = Nu_e_vec(71);

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

%%%%% low-hybrid frequency range
N = 80;        %%% must be even
N_upper = 400; %%% must be even
upper_Bound = 32000;
q_n0 = [2; 20; 800; 16000];


%%%%% upper-hybrid frequency range
N = 100;       %%% must be even
N_upper = 400; %%% must be even
upper_Bound = 41;
q_n0 = [2; 4];
[dq_simp, q_cs, p_cs] = quadratureMethod_forIntegralEqs_multipleIntervals('difficult_simpsonMarkovGildendurgExpWithCollisions', N, N_upper, q_n0, 1e-8, upper_Bound);



upper_Bound = 30;
N = 20;         %%% must be even
N_upper = 400;  %%% must be even
q_n00 = sort(real(q1(1:30)));
q_n00 = [1.035; 1.134; q_n00(1); q_n00(3); q_n00(5); q_n00(7:end)];
stepQ = 0.01;
q_n0 = [];
for iq0 = 1:5
    q_n0 = [q_n0; q_n00(iq0)-stepQ; q_n00(iq0); q_n00(iq0)+stepQ];
end
stepQ = 0.1;
for iq0 = 6:size(q_n00,1)
    q_n0 = [q_n0; q_n00(iq0)-stepQ; q_n00(iq0); q_n00(iq0)+stepQ];
end
[dq_simp, q_cs, p_cs] = quadratureMethod_forIntegralEqs_multipleIntervals('difficult_simpsonMarkovGildendurgExpWithCollisions', N, N_upper, q_n0, 1e-8, upper_Bound);




% upper_Bound = 35;
% N = 1600;
% [dq_simp, q_cs, p_cs] = quadratureMethod_forIntegralEqs('simpson', N, 1i, 1e-8, upper_Bound);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%% test psy %%%
% psi = psi1_q__gyrotropic(k_0, k_0, a_0, EE1, GG1, HH1, m, p_cs, q_cs);
% psi2 = psi2_q__gyrotropic(k_0, k_0, a_0, EE1, GG1, HH1, m, p_cs, q_cs);
% figure(112);
% hold on
% plot(q_cs, real(psi))
% plot(q_cs, real(psi2))
% hold off
% %%%%%%%%%%%%%%%%%


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

interalQ = q_cs >= 1;
% a_cs_alfa1_forw(interalQ) = 0*a_cs_alfa1_forw(interalQ);
% a_cs_alfa2_forw(interalQ) = 0*a_cs_alfa2_forw(interalQ);

% a_cs_alfa1_forw(isinf(a_cs_alfa1_forw)) = isnan(a_cs_alfa1_forw(isinf(a_cs_alfa1_forw))) * 0;
% a_cs_alfa2_forw(isinf(a_cs_alfa2_forw)) = isnan(a_cs_alfa2_forw(isinf(a_cs_alfa2_forw))) * 0;
a_cs_alfa1_forw(isinf(a_cs_alfa1_forw)) = 0;
a_cs_alfa2_forw(isinf(a_cs_alfa2_forw)) = 0;
a_cs_alfa1_forw(isnan(a_cs_alfa1_forw)) = 0;
a_cs_alfa2_forw(isnan(a_cs_alfa2_forw)) = 0;

% a_cs_alfa1_forw(interalQ) = 0*a_cs_alfa1_forw(interalQ);
% a_cs_alfa2_forw(interalQ) = 0*a_cs_alfa2_forw(interalQ);

% a_cs_alfa1_forw = 0*a_cs_alfa1_forw;
% a_cs_alfa2_forw = 0*a_cs_alfa2_forw;
% a_smn_forw = a_smn_forw*0;
% a_cs_alfa1_forw(end-300) = 1;
% a_cs_alfa2_forw(end-300) = 1;



% figure(113)
% plot(q_cs,abs(a_cs_alfa1_forw));
% figure(114)
% plot(q_cs,abs(a_cs_alfa2_forw));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a_cs_alfa1_back(isinf(a_cs_alfa1_back)) = isnan(a_cs_alfa1_back(isinf(a_cs_alfa1_back))) * 0;
% a_cs_alfa2_back(isinf(a_cs_alfa2_back)) = isnan(a_cs_alfa2_back(isinf(a_cs_alfa2_back))) * 0;
a_cs_alfa1_back(isinf(a_cs_alfa1_back)) = 0;
a_cs_alfa2_back(isinf(a_cs_alfa2_back)) = 0;
a_cs_alfa1_back(isnan(a_cs_alfa1_back)) = 0;
a_cs_alfa2_back(isnan(a_cs_alfa2_back)) = 0;

a_cs_alfa1_back_prop = a_cs_alfa1_back;
a_cs_alfa2_back_prop = a_cs_alfa2_back;
a_cs_alfa1_back_prop(interalQ) = 0*a_cs_alfa1_back_prop(interalQ);
a_cs_alfa2_back_prop(interalQ) = 0*a_cs_alfa2_back_prop(interalQ);

a_cs_alfa1_forw_prop = a_cs_alfa1_forw;
a_cs_alfa2_forw_prop = a_cs_alfa2_forw;
a_cs_alfa1_forw_prop(interalQ) = 0*a_cs_alfa1_forw_prop(interalQ);
a_cs_alfa2_forw_prop(interalQ) = 0*a_cs_alfa2_forw_prop(interalQ);

z =  2*pi / k_0 * (0);
k_0 = sourceParameters.k_0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% forward excited waves %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%L = sourceParameters.zCoordinate;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda_0 = 2 * pi / k_0;
L_0 = 4.1*a_0; L_MID = 0.1*lambda_0; L_END = 0.5*lambda_0;
NL1 = 500; % % требуетс€ как минимум 1120 точек на одну длину волны
NL2 = 500;
LL1 = L_0:(L_MID - L_0)/NL1:L_MID;
LL2 = L_MID:(L_END - L_MID)/NL2:L_END;
% LL = [LL1 LL2];
LL = 4*d;%0.5*lambda_0;%[0.2 0.5 1.5]*lambda_0;
LL = 4 * d;
P_mod_back = zeros(size(LL,2),1);
P_cs_back  = zeros(size(LL,2),1);
P_cs_forw  = zeros(size(LL,2),1);
for il = 1:size(LL,2) % new cycle on L created
    %%% only for upper-hybrid range! %%%%%%%%%%
%     clear p_n
%     p_n__of_descreteMode_of_gyrotropicCyl
%     q_n = sqrt(1-p_n.^2);
%     q_n = q_n.* (2*(imag(q_n) <= 0)-1);
    
    
    L = - LL(il);
    sourceParameters.zCoordinate = - LL(il);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



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

                      
figure(113)
hold on;
     plot(q_cs,abs(b_p_field_1_back));
     plot(real(q1),0*q1, 'o');
     plot(real(q2),0*q2, '*');
hold off;
figure(114)
hold on;
     plot(q_cs,abs(b_p_field_2_back));
     plot(real(q1),0*q1, 'o');
     plot(real(q2),0*q2, '*');
hold off;


%%%%%%%% CALCULATING OF PARTIAL POWERS %%%%%%%%%%% added by Oleg
%%%%%%%% Ostafiychuk 14.01.2020
%clc

b_p_field_1_back_prop = [b_p_field_1_back(~interalQ) 0*b_p_field_1_back(interalQ)];
b_p_field_2_back_prop = [b_p_field_2_back(~interalQ) 0*b_p_field_2_back(interalQ)];

a_p_Ebeam_forw_prop = [a_p_Ebeam_forw(~interalQ) 0*a_p_Ebeam_forw(interalQ)];
a_p_Hbeam_forw_prop = [a_p_Hbeam_forw(~interalQ) 0*a_p_Hbeam_forw(interalQ)];

% this is only for UH range !!! %%%%%%%%%%%%%%%%%%%%%%
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
%     P_n_minus(in) = P_of_discreteMode(typeOfCylinder, q_n(in), - p_n(in), k_0, a_0, ...
%                                                 EE1, GG1, HH1, MU1, EE, MU, m, a_smn_back(in).*exp(1i * k_0 * p_n(in) * abs(L)) + b_n_back(in));

      P_n_minus(in) = P_of_discreteMode_NEW( - p_n(in), waveguideParameters, sourceParameters,...
          a_smn_back(in).*exp(1i * k_0 *  p_n(in) * abs(L)) + b_n_back(in));
    
end

P_mod_back(il) = - sum(P_n_minus);

Phase_exp_cs = exp(1i * k_0 * p_cs * abs(L));
%Phase_exp_cs_beta = exp(1i * k_0 * p_cs * abs(L - d));
Phase_exp_cs(interalQ) = 0;
%Phase_exp_cs_beta(interalQ) = 0;
% P_cs_minus_1 = ...
%       sum(dq_simp.*P_of_continuousWaves_alpha1(typeOfCylinder, q_cs, p_0, - p_cs, k_0, a_0, ...
%                                  EE1, GG1, HH1, MU1, EE, MU, c, m, z, j_f, j_z, d, a_cs_alfa1_back_prop.* Phase_exp_cs  + b_p_field_1_back_prop));
% P_cs_minus_2 = ...
%       sum(dq_simp.*P_of_continuousWaves_alpha2(typeOfCylinder, q_cs, p_0, - p_cs, k_0, a_0, ...
%                                  EE1, GG1, HH1, MU1, EE, MU, c, m, z, j_f, j_z, d, a_cs_alfa2_back_prop.* Phase_exp_cs + b_p_field_2_back_prop));
P_cs_minus_1 = ...
      sum(dq_simp.*P_of_continuousWaves_alpha1(typeOfCylinder, q_cs, p_0, - p_cs, k_0, a_0, ...
                                 EE1, GG1, HH1, MU1, EE, MU, c, m, z, j_f, j_z, d,...
                                 a_cs_alfa1_back_prop.* Phase_exp_cs  + b_p_field_1_back_prop));
P_cs_minus_2 = ...
      sum(dq_simp.*P_of_continuousWaves_alpha2(typeOfCylinder, q_cs, p_0, - p_cs, k_0, a_0, ...
                                 EE1, GG1, HH1, MU1, EE, MU, c, m, z, j_f, j_z, d,...
                                 a_cs_alfa2_back_prop.* Phase_exp_cs + b_p_field_2_back_prop));

P_cs_back(il) = - P_cs_minus_1 - P_cs_minus_2;

%%%%%%%% in half-space z>0%%%%%%%%%%

%%% the norms of eigenwaves of free space %%%
% N_p1 = @(q, p) (-1)^(m+1) * c * p./ (k_0.^2 * q);
% N_p2 = @(q, p) (-1)^(m+2) * c * p./ (k_0.^2 * q);
% P_cs_plus_1 = 1/4.*sum(dq_simp.*abs(a_p_Ebeam_forw_prop).^2.*N_p1(q_cs, p_cs));
% P_cs_plus_2 = 1/4.*sum(dq_simp.*abs(a_p_Hbeam_forw_prop).^2.*N_p2(q_cs, p_cs));

%P_p1 = @(q, p) 1/4.*(-1)^(m+1) * c * p./ (k_0.^2 * q); %%%%% <- probably, incorrect formulas
%P_p2 = @(q, p) 1/4.*(-1)^(m+2) * c * p./ (k_0.^2 * q)*(-1);

P_p1 = @(q, p) 1/4. * c * p./ (k_0.^2 * q);
P_p2 = @(q, p) 1/4. * c * p./ (k_0.^2 * q);

P_cs_plus_1 = sum(dq_simp.*abs(a_p_Ebeam_forw_prop).^2.*P_p1(q_cs, p_cs));
P_cs_plus_2 = sum(dq_simp.*abs(a_p_Hbeam_forw_prop).^2.*P_p2(q_cs, p_cs));
                             
P_cs_forw(il) = P_cs_plus_1 + P_cs_plus_2;   

P_ratio = P_mod_back/(P_cs_back + P_cs_forw);

% figure(10)
% bar(real(p_n), - P_n_minus);
% Ferr = getframe;

end

%% %%% in the presence of an infinite waveguide <-- new added! 
P_n_minus_waveguide = zeros(size(q_n,1),1);
P_n_plus_waveguide = zeros(size(q_n,1),1);
% sourceParameters.zCoordinate = 0;
for in = 1:size(q_n,1)
%     P_n_minus_waveguide(in) = P_of_discreteMode_NEW(typeOfCylinder, q_n(in), - p_n(in), k_0, a_0, ...
%                                                 EE1, GG1, HH1, MU1, EE, MU, m, a_smn_back(in));
%     
%     P_n_plus_waveguide(in) = P_of_discreteMode_NEW(typeOfCylinder, q_n(in), p_n(in), k_0, a_0, ...
%                                                 EE1, GG1, HH1, MU1, EE, MU, m, a_smn_forw(in));   
    
    P_n_minus_waveguide(in) = P_of_discreteMode_NEW( - p_n(in), waveguideParameters, sourceParameters, a_smn_back(in));
    
    P_n_plus_waveguide(in) = P_of_discreteMode_NEW( p_n(in), waveguideParameters, sourceParameters, a_smn_forw(in));   
end

P_mod_back_waveguide = - sum(P_n_minus_waveguide);
P_mod_forw_waveguide =   sum(P_n_plus_waveguide);
P_mod_waveguide = P_mod_back_waveguide + P_mod_forw_waveguide;
% bar(-p_n,P_n_plus_waveguide-P_n_minus_waveguide);

P_cs_minus_1_waveguide = ...
      sum(dq_simp.*P_of_continuousWaves_alpha1(typeOfCylinder, q_cs, p_0, - p_cs, k_0, a_0, ...
                                 EE1, GG1, HH1, MU1, EE, MU, c, m, z, j_f, j_z, d, a_cs_alfa1_back_prop));
P_cs_minus_2_waveguide = ...
      sum(dq_simp.*P_of_continuousWaves_alpha2(typeOfCylinder, q_cs, p_0, - p_cs, k_0, a_0, ...
                                 EE1, GG1, HH1, MU1, EE, MU, c, m, z, j_f, j_z, d, a_cs_alfa2_back_prop));
P_cs_plus_1_waveguide = ...
      sum(dq_simp.*P_of_continuousWaves_alpha1(typeOfCylinder, q_cs, p_0, p_cs, k_0, a_0, ...
                                 EE1, GG1, HH1, MU1, EE, MU, c, m, z, j_f, j_z, d, a_cs_alfa1_forw_prop));   
P_cs_plus_2_waveguide = ...
      sum(dq_simp.*P_of_continuousWaves_alpha2(typeOfCylinder, q_cs, p_0, p_cs, k_0, a_0, ...
                                 EE1, GG1, HH1, MU1, EE, MU, c, m, z, j_f, j_z, d, a_cs_alfa2_forw_prop));                             

% pp = @(q) sqrt(1 - q.^2);
% P_cs_minus_1_waveguide = quadgk(@(q)P_of_continuousWaves_alpha1(typeOfCylinder, q, p_0, - pp(q), k_0, a_0, ...
%                                  EE1, GG1, HH1, MU1, EE, MU, c, m, z, j_f, j_z, d, ...
%                                  a_sma_of_continuousWaves_alpha1_Decorator(typeOfCylinder, q,  - pp(q), waveguideParameters, sourceParameters)),0,1);
%                              
% P_cs_minus_2_waveguide = quadgk(@(q)P_of_continuousWaves_alpha2(typeOfCylinder, q, p_0, - pp(q), k_0, a_0, ...
%                                  EE1, GG1, HH1, MU1, EE, MU, c, m, z, j_f, j_z, d, ...
%                                  a_sma_of_continuousWaves_alpha2_Decorator(typeOfCylinder, q,  - pp(q), waveguideParameters, sourceParameters)),0,1);
%                              
% P_cs_plus_1_waveguide = quadgk(@(q)P_of_continuousWaves_alpha1(typeOfCylinder, q, p_0, pp(q), k_0, a_0, ...
%                                  EE1, GG1, HH1, MU1, EE, MU, c, m, z, j_f, j_z, d, ...
%                                  a_sma_of_continuousWaves_alpha1_Decorator(typeOfCylinder, q,  pp(q), waveguideParameters, sourceParameters)),0,1);                             
%                             
% P_cs_plus_2_waveguide = quadgk(@(q)P_of_continuousWaves_alpha2(typeOfCylinder, q, p_0, pp(q), k_0, a_0, ...
%                                  EE1, GG1, HH1, MU1, EE, MU, c, m, z, j_f, j_z, d, ...
%                                  a_sma_of_continuousWaves_alpha2_Decorator(typeOfCylinder, q,  pp(q), waveguideParameters, sourceParameters)),0,1);
                              
P_cs_back_waveguide = - P_cs_minus_1_waveguide - P_cs_minus_2_waveguide;
P_cs_forw_waveguide =   P_cs_plus_1_waveguide + P_cs_plus_2_waveguide;

P_cs_waveguide = P_cs_back_waveguide + P_cs_forw_waveguide;
R_cs_waveguide = P_cs_waveguide*2 / I_0^2 * toOhms

R_mod_waveguide = P_mod_waveguide*2 / I_0^2 * toOhms

P_sum_waveguide = R_mod_waveguide + P_cs_waveguide;

%% calculation of the power radiated by the source located in free space
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
P_free_space = P_free_space_E + P_free_space_H;
R_free_space = P_free_space*2 / I_0^2 * toOhms;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%P_sum_waveguide = 32.9754;
figure(7)
plot(LL/lambda_0, P_mod_back/P_sum_waveguide, 'r'); hold on
plot(LL/lambda_0, 10*P_cs_back/P_sum_waveguide, 'b'); hold on
plot(LL/lambda_0, 10*P_cs_forw/P_sum_waveguide, 'k'); hold on

%P_Ratio = P_mod_back/(P_cs_back + P_cs_forw)

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

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
Ez_forwardExcitedWave = zeros(1,size(rho,2));
Ephi_forwardExcitedWave = zeros(1,size(rho,2));
Hphi_forwardExcitedWave = zeros(1,size(rho,2));
%%% discrete mode
for in = 1:size(q_n,1)
    Ez_forwardExcitedWave = Ez_forwardExcitedWave +...
        field_of_discreteMode_Decorator('Ez',   typeOfCylinder, rho, z, q_n(in), p_n(in), a_smn_forw(in), waveguideParameters, sourceParameters).*...
        exp(1i * k_0 * p_n(in) * L);

    Ephi_forwardExcitedWave = Ephi_forwardExcitedWave +...
        field_of_discreteMode_Decorator('Ephi',   typeOfCylinder, rho, z, q_n(in), p_n(in), a_smn_forw(in), waveguideParameters, sourceParameters).*...
        exp(1i * k_0 * p_n(in) * L);
    
    Hphi_forwardExcitedWave = Hphi_forwardExcitedWave +...
        field_of_discreteMode_Decorator('Hphi',   typeOfCylinder, rho, z, q_n(in), p_n(in), a_smn_forw(in), waveguideParameters, sourceParameters).*...
        exp(1i * k_0 * p_n(in) * L);
end


%%% continues spectrume
Expk0pcsL = exp(1i * k_0 * p_cs * L);
a_cs_alfa1_forw_Exp = a_cs_alfa1_forw.* exp(1i * k_0 * p_cs * L);
a_cs_alfa2_forw_Exp = a_cs_alfa2_forw.* exp(1i * k_0 * p_cs * L);
for ir = 1:size(rho,2)
    Ez_forwardExcitedWave(ir) = Ez_forwardExcitedWave(ir) +...
          (k_0/2/pi) * sum(dq_simp.* field__ContinuesWaves_discreteRepr('Ez', typeOfCylinder, rho(ir), z, q_cs, p_cs,...
          a_cs_alfa1_forw_Exp, a_cs_alfa1_forw_Exp*0, a_cs_alfa2_forw_Exp, a_cs_alfa2_forw_Exp*0,...
          waveguideParameters, sourceParameters, coeffsOfField));
      
    Ephi_forwardExcitedWave(ir) = Ephi_forwardExcitedWave(ir) +...
          (k_0/2/pi) * sum(dq_simp.* field__ContinuesWaves_discreteRepr('Ephi', typeOfCylinder, rho(ir), z, q_cs, p_cs,...
          a_cs_alfa1_forw_Exp, a_cs_alfa1_forw_Exp*0, a_cs_alfa2_forw_Exp, a_cs_alfa2_forw_Exp*0,...
          waveguideParameters, sourceParameters, coeffsOfField));
      
    Hphi_forwardExcitedWave(ir) = Hphi_forwardExcitedWave(ir) +...
          (k_0/2/pi) * sum(dq_simp.* field__ContinuesWaves_discreteRepr('Hphi', typeOfCylinder, rho(ir), z, q_cs, p_cs,...
          a_cs_alfa1_forw_Exp, a_cs_alfa1_forw_Exp*0, a_cs_alfa2_forw_Exp, a_cs_alfa2_forw_Exp*0,...
          waveguideParameters, sourceParameters, coeffsOfField));
end

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

Ez_waveguidespace = Ez_forwardExcitedWave;
Ephi_waveguidespace = Ephi_forwardExcitedWave;
Hphi_waveguidespace = Hphi_forwardExcitedWave;

%%% add field field of diffracted waves - discrete spectrum waves
for in = 1:size(q_n,1)
    Ez_waveguidespace = Ez_waveguidespace +...
        field_of_discreteMode_Decorator('Ez',   typeOfCylinder, rho, z, q_n(in), -p_n(in), b_n_back(in), waveguideParameters, sourceParameters);

    Ephi_waveguidespace = Ephi_waveguidespace +...
        field_of_discreteMode_Decorator('Ephi',   typeOfCylinder, rho, z, q_n(in), -p_n(in), b_n_back(in), waveguideParameters, sourceParameters);
    
    Hphi_waveguidespace = Hphi_waveguidespace +...
        field_of_discreteMode_Decorator('Hphi',   typeOfCylinder, rho, z, q_n(in), -p_n(in), b_n_back(in), waveguideParameters, sourceParameters);
end   

%%% add field field of diffracted waves - continuous spectrum waves
for ir = 1:size(rho,2)
    Ez_waveguidespace(ir) = Ez_waveguidespace(ir) +...
          (k_0/2/pi) * sum(dq_simp.* field__ContinuesWaves_discreteRepr('Ez', typeOfCylinder, rho(ir), z, q_cs, p_cs,...
          b_p_field_1_back*0, b_p_field_1_back, b_p_field_2_back*0, b_p_field_2_back,...
          waveguideParameters, sourceParameters, coeffsOfField));
    
    Ephi_waveguidespace(ir) = Ephi_waveguidespace(ir) +...
          (k_0/2/pi) * sum(dq_simp.* field__ContinuesWaves_discreteRepr('Ephi', typeOfCylinder, rho(ir), z, q_cs, p_cs,...
          b_p_field_1_back*0, b_p_field_1_back, b_p_field_2_back*0, b_p_field_2_back,...
          waveguideParameters, sourceParameters, coeffsOfField));
     
    Hphi_waveguidespace(ir) = Hphi_waveguidespace(ir) +...
          (k_0/2/pi) * sum(dq_simp.* field__ContinuesWaves_discreteRepr('Hphi', typeOfCylinder, rho(ir), z, q_cs, p_cs,...
          b_p_field_1_back*0, b_p_field_1_back, b_p_field_2_back*0, b_p_field_2_back,...
          waveguideParameters, sourceParameters, coeffsOfField));
end

figure(1)
hold on
plot(rho / a_0, real(Ez_waveguidespace), 'b-o', rho / a_0, imag(Ez_waveguidespace),'r-o')
title('E_{z}');
hold off
figure(2)
hold on
plot(rho / a_0, real(Ephi_waveguidespace), 'b-o', rho / a_0, imag(Ephi_waveguidespace),'r-o')
title('E_{\phi}');
hold off
figure(3)
hold on
plot(rho / a_0, real(Hphi_waveguidespace), 'b-o', rho / a_0, imag(Hphi_waveguidespace),'r-o')
title('H_{\phi}');
hold off

% %%%%%%%%%%% For free space part %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_p_Ebeam_freespace = a_p_Ebeam_forw;
a_p_Hbeam_freespace = a_p_Hbeam_forw;
    
Ez_freespace   = zeros(size(rho));
Ephi_freespace = zeros(size(rho));
% Erho_freespace = zeros(size(rho));
% Hz_freespace   = zeros(size(rho));
Hphi_freespace = zeros(size(rho));
% Hrho_freespace = zeros(size(rho));
for ir = 1:size(rho,2)
    Ez_freespace(ir) = Ez_freespace(ir) +...
          (k_0/2/pi) * sum(dq_simp.* field__ForBeamDiffraction_InFreeSpace_discreteRepr('Ez', typeOfCylinder,rho(ir), z, q_cs, p_cs, k_0,  m,...
          a_p_Ebeam_freespace, a_p_Hbeam_freespace));
      
    Ephi_freespace(ir) = Ephi_freespace(ir) +...
          (k_0/2/pi) * sum(dq_simp.* field__ForBeamDiffraction_InFreeSpace_discreteRepr('Ephi', typeOfCylinder,rho(ir), z, q_cs, p_cs, k_0,  m,...
          a_p_Ebeam_freespace, a_p_Hbeam_freespace));
              
    Hphi_freespace(ir) = Hphi_freespace(ir) +...
          (k_0/2/pi) * sum(dq_simp.* field__ForBeamDiffraction_InFreeSpace_discreteRepr('Hphi', typeOfCylinder,rho(ir), z, q_cs, p_cs, k_0,  m,...
          a_p_Ebeam_freespace, a_p_Hbeam_freespace));
end

figure(1)
hold on
plot(rho / a_0, real(Ez_freespace), 'b', rho / a_0, imag(Ez_freespace),'r')
title('E_{z}');
hold off
figure(2)
hold on
plot(rho / a_0, real(Ephi_freespace), 'b',rho / a_0, imag(Ephi_freespace),'r')
title('E_{\phi}');
hold off
figure(3)
hold on
plot(rho / a_0, real(Hphi_freespace), 'b',rho / a_0, imag(Hphi_freespace),'r')
title('H_{\phi}');
hold off

toc
  
  
