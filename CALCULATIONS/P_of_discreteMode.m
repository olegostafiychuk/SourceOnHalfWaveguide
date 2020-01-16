function P = P_of_discreteMode(typeOfCylinder, q, p, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, aplus2)
   %%%%%% norm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    N_n = Norm_of_descreteMode(typeOfCylinder, q, p, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m);
    
   %%%%%% power %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    P = 1/4 * real(abs(aplus2.*conj(aplus2)).* N_n); 



    
    
    
    
    
    
    
    
    
