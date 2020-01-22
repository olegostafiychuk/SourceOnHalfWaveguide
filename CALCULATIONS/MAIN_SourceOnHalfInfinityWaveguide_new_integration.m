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
end
q_n = sqrt(1-p_n.^2);
q_n = q_n.* (2*(imag(q_n) <= 0)-1);

%upper_Bound = 32000;
%upper_Bound = 16000;
% upper_Bound = 8000;
%N = 1600;
% N = 3200;

% upper_Bound = 400;
% N = 40;

%[dq_simp, q_cs, p_cs] = quadratureMethod_forIntegralEqs('simpson', N, 1i, 1e-8, upper_Bound);
% [dq_simp, q_cs, p_cs] = quadratureMethod_forIntegralEqs('difficult_simpson', N, q_0, 1e-7, upper_Bound);



N = 40;         %%% must be even
N_upper = 800; %%% must be even
upper_Bound = 24000;
q_n0 = [2; 20; 400; upper_Bound];
[dq_simp, q_cs, p_cs] = quadratureMethod_forIntegralEqs_multipleIntervals('difficult_simpsonMarkovGildendurgExpWithCollisions', N, N_upper, q_n0, 1e-8, upper_Bound);




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
% interalQ = q_cs > 0.99;
% a_cs_alfa1_forw(interalQ) = 0*a_cs_alfa1_forw(interalQ);
% a_cs_alfa2_forw(interalQ) = 0*a_cs_alfa2_forw(interalQ);

% a_cs_alfa1_forw(1) = 1;
% a_cs_alfa1_forw(isnan(a_cs_alfa1_forw)) = isnan(a_cs_alfa1_forw(isnan(a_cs_alfa1_forw))) * 0;
% a_cs_alfa2_forw(isnan(a_cs_alfa2_forw)) = isnan(a_cs_alfa2_forw(isnan(a_cs_alfa2_forw))) * 0;

a_cs_alfa1_forw(isinf(a_cs_alfa1_forw)) = isnan(a_cs_alfa1_forw(isinf(a_cs_alfa1_forw))) * 0;
a_cs_alfa2_forw(isinf(a_cs_alfa2_forw)) = isnan(a_cs_alfa2_forw(isinf(a_cs_alfa2_forw))) * 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% forward excited waves %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z =  2*pi / k_0 * (0);
L = sourceParameters.zCoordinate;
k_0 = sourceParameters.k_0;
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
  
  
