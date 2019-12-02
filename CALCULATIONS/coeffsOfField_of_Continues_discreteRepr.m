%%%%% расчёт коэффициента для волн непрерывного спектра в прямом
%%%%% направлении (в сторону среды 2)
%%%%% второе приближение в нахождении коэффициента aplus
function coeffsOfField =...
          coeffsOfField_of_Continues_discreteRepr(typeOfCylinder, q, p, waveguideParameters, sourceParameters)

m = sourceParameters.m;
k_0 = sourceParameters.k_0;
k = k_0;

EE1 = waveguideParameters.EEinner;
GG1 = waveguideParameters.GGinner;
HH1 = waveguideParameters.HHinner;
MU1 = waveguideParameters.MUinner;
EE = waveguideParameters.EEouter;
MU = waveguideParameters.MUouter;
a_0 = waveguideParameters.radius;
      
%   p = sqrt(1-q.^2);
%   p = real(p) - 1i * abs(imag(p));
    
  m0 = m;
  
switch(typeOfCylinder)
    case 'PerfectConductivity'
        
        m = -m;
        
        %%%%%%%% forward wave %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%% first type of waves %%%%%%%%%%%%%%%%%%%%%%%%%%%        
        psi1 = -besselh(m,2, k_0 * q * a_0)./ besselh(m,1, k_0 * q * a_0);
        Dm2_forward_1 = 0;
        Cm2_forward_1 = 1;
    
        
        
        
        %%%%%%%%%%%%%%%%%% second type of waves %%%%%%%%%%%%%%%%%%%%%%%%%%%
        psi2 =  (besselh(m,2, k_0 * q * a_0) * m./ (k_0 * q * a_0) - besselh(m+1,2, k_0 * q * a_0))./...
                (besselh(m,1, k_0 * q * a_0) * m./ (k_0 * q * a_0) - besselh(m+1,1, k_0 * q * a_0));
        Dm2_forward_2 = 1;
        Cm2_forward_2 = 0;
        
        psi_forward_1  = psi1;
        psi_backward_1 = psi1;
        psi_forward_2  = psi2;
        psi_backward_2 = psi2;
        psi_forward_1_transp  = [];
        psi_backward_1_transp = [];
        psi_forward_2_transp  = [];
        psi_backward_2_transp = [];
        
        B_1_forward_1 = 0;B_1_forward_2 = 0;
        B_2_forward_1 = 0;B_2_forward_2 = 0;
        B_1_backward_1 = 0;B_1_backward_2 = 0;
        B_2_backward_1 = 0;B_2_backward_2 = 0;
        Cm2_backward_1 = Cm2_forward_1;
        Cm2_backward_2 = Cm2_forward_2;
        Dm2_backward_1 = Dm2_forward_1;
        Dm2_backward_2 = Dm2_forward_2;
        
        
    case 'Isotropic'      
        %%% вычисляем первый тип волн непрерывного спектра
       
        %%% first type of waves
        psi_forward_1  = psi1_q(k_0, k, a_0, EE1, MU1, EE, MU, m,   p, q);
        psi_backward_1  = psi1_q(k_0, k, a_0, EE1, MU1, EE, MU, m,  -p, q);
        [B_1_forward_1,B_2_forward_1,Cm2_forward_1,  Dm2_forward_1]  = coefficientsOfContinuousSpectrum(typeOfCylinder,k_0, k, a_0, EE1,0,0, MU1, EE, MU, m,   p, q,   psi_forward_1,  1);
        [B_1_backward_1,B_2_backward_1,Cm2_backward_1, Dm2_backward_1] = coefficientsOfContinuousSpectrum(typeOfCylinder,k_0, k, a_0, EE1,0,0, MU1, EE, MU, m, (-p), q, psi_backward_1, 1);
%         [B_1_back_PandM_1,B_2_back_PandM_1,Cm2_back_PandM_1, Dm2_back_PandM_1] = coefficientsOfContinuousSpectrum(typeOfCylinder,k_0, k, a_0, EE1, MU1, EE, MU, -m, (-p), q, psi_forward_1, 1);
        
        %%% second type of waves
        psi_forward_2  = psi2_q(k_0, k, a_0, EE1, MU1, EE, MU, m,   p, q);
        psi_backward_2  = psi2_q(k_0, k, a_0, EE1, MU1, EE, MU, m,  -p, q);
        [B_1_forward_2,B_2_forward_2,Cm2_forward_2,  Dm2_forward_2]  = coefficientsOfContinuousSpectrum_forPsi2(typeOfCylinder,k_0, k, a_0, EE1,0,0, MU1, EE, MU, m,   p, q,   psi_forward_2,  1);
        [B_1_backward_2,B_2_backward_2,Cm2_backward_2, Dm2_backward_2] = coefficientsOfContinuousSpectrum_forPsi2(typeOfCylinder,k_0, k, a_0, EE1,0,0, MU1, EE, MU, m, (-p), q, psi_backward_2, 1);
%         [B_1_back_PandM_2,B_2_back_PandM_2,Cm2_back_PandM_2, Dm2_back_PandM_2] = coefficientsOfContinuousSpectrum_forPsi2(typeOfCylinder,k_0, k, a_0, EE1, MU1, EE, MU, -m, (-p), q, psi_forward_2, 1);
        
        psi_forward_1_transp  = [];
        psi_backward_1_transp = [];
        psi_forward_2_transp  = [];
        psi_backward_2_transp = [];
        
    case 'Gyrotropic'      
        %%% вычисляем первый тип волн непрерывного спектра
       
        %%% first type of waves
        psi_forward_1   = psi1_q__gyrotropic(k_0, k, a_0, EE1, GG1, HH1, m,   p, q);
        psi_backward_1  = psi1_q__gyrotropic(k_0, k, a_0, EE1, GG1, HH1, m,  -p, q);
        psi_forward_1_transp   = psi1_q__gyrotropic(k_0, k, a_0, EE1, -GG1, HH1, -m,   p, q);
        psi_backward_1_transp  = psi1_q__gyrotropic(k_0, k, a_0, EE1, -GG1, HH1, -m,  -p, q);
%         psi_forward_1_transp   = psi_backward_1;
%         psi_backward_1_transp  = psi_forward_1;
        [B_1_forward_1,B_2_forward_1,Cm2_forward_1,  Dm2_forward_1]  = coefficientsOfContinuousSpectrum(typeOfCylinder,k_0, k, a_0, EE1, GG1, HH1, MU1, EE, MU, m,   p, q,   psi_forward_1, 1);
        [B_1_backward_1,B_2_backward_1,Cm2_backward_1, Dm2_backward_1] = coefficientsOfContinuousSpectrum(typeOfCylinder,k_0, k, a_0, EE1, GG1, HH1, MU1, EE, MU, m, (-p), q, psi_backward_1, 1);
%         [B_1_back_PandM_1,B_2_back_PandM_1,Cm2_back_PandM_1, Dm2_back_PandM_1] = coefficientsOfContinuousSpectrum(typeOfCylinder,k_0, k, a_0, EE1, MU1, EE, MU, -m, (-p), q, psi_forward_1, 1);
        
        %%% second type of waves
        psi_forward_2  = psi2_q__gyrotropic(k_0, k, a_0, EE1, GG1, HH1, m,   p, q);
        psi_backward_2  = psi2_q__gyrotropic(k_0, k, a_0, EE1, GG1, HH1, m,  -p, q);
        psi_forward_2_transp   = psi2_q__gyrotropic(k_0, k, a_0, EE1, -GG1, HH1, -m,   p, q);
        psi_backward_2_transp  = psi2_q__gyrotropic(k_0, k, a_0, EE1, -GG1, HH1, -m,  -p, q);
%         psi_forward_2_transp   = psi_backward_2;
%         psi_backward_2_transp  = psi_forward_2;
        [B_1_forward_2,B_2_forward_2,Cm2_forward_2,  Dm2_forward_2]  = coefficientsOfContinuousSpectrum_forPsi2(typeOfCylinder,k_0, k, a_0, EE1, GG1, HH1, MU1, EE, MU, m,   p, q,  psi_forward_2,  1);
        [B_1_backward_2,B_2_backward_2,Cm2_backward_2, Dm2_backward_2] = coefficientsOfContinuousSpectrum_forPsi2(typeOfCylinder,k_0, k, a_0, EE1, GG1, HH1, MU1, EE, MU, m, (-p), q, psi_backward_2, 1);
%         [B_1_back_PandM_2,B_2_back_PandM_2,Cm2_back_PandM_2, Dm2_back_PandM_2] = coefficientsOfContinuousSpectrum_forPsi2(typeOfCylinder,k_0, k, a_0, EE1, MU1, EE, MU, -m, (-p), q, psi_forward_2, 1);

end

coeffsOfField.psi_forward_1 = psi_forward_1;
coeffsOfField.psi_backward_1 = psi_backward_1;
coeffsOfField.psi_forward_1_transp = psi_forward_1_transp;
coeffsOfField.psi_backward_1_transp = psi_backward_1_transp;
coeffsOfField.psi_forward_2 = psi_forward_2;
coeffsOfField.psi_backward_2 = psi_backward_2;
coeffsOfField.psi_forward_2_transp = psi_forward_2_transp;
coeffsOfField.psi_backward_2_transp = psi_backward_2_transp;
coeffsOfField.B_1_forward_1 = B_1_forward_1;
coeffsOfField.B_2_forward_1 = B_2_forward_1;
coeffsOfField.Cm2_forward_1 = Cm2_forward_1;
coeffsOfField.Dm2_forward_1 = Dm2_forward_1;
coeffsOfField.B_1_backward_1 = B_1_backward_1;
coeffsOfField.B_2_backward_1 = B_2_backward_1;
coeffsOfField.Cm2_backward_1 = Cm2_backward_1;
coeffsOfField.Dm2_backward_1 = Dm2_backward_1;
coeffsOfField.B_1_forward_2 = B_1_forward_2;
coeffsOfField.B_2_forward_2 = B_2_forward_2;
coeffsOfField.Cm2_forward_2 = Cm2_forward_2;
coeffsOfField.Dm2_forward_2 = Dm2_forward_2;
coeffsOfField.B_1_backward_2 = B_1_backward_2;
coeffsOfField.B_2_backward_2 = B_2_backward_2;
coeffsOfField.Cm2_backward_2 = Cm2_backward_2;
coeffsOfField.Dm2_backward_2 = Dm2_backward_2;
