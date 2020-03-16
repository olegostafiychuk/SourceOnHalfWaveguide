function P = P_of_discreteMode_NEW(p, waveguideParameters, sourceParameters, a_mn)
   %%%%%% norm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     N_n = Norm_of_descreteMode(typeOfCylinder, q, p, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m);
    [Pn, innerAnalytical, outerAnalytical] = powerOfDiscreteSpectrumMode(p, waveguideParameters, sourceParameters);
    
   %%%%%% power %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     P = 1/4 * real(abs(aplus2.*conj(aplus2)).* N_n); 
P = real(abs(a_mn.*conj(a_mn)).* Pn); 
