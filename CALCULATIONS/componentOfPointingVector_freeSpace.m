function S_component = componentOfPointingVector_freeSpace(componentOfPV, typeOfCylinder, rho, phi, z, q, p, a_excitationOfEwave, a_excitationOfHwave,...
            k_0, m, c, dq_simp)
        
        switch(componentOfPV)
             case 'S_rho'
                    E_phi = (k_0/2/pi) * sum(dq_simp.* field__ForBeamDiffraction_InFreeSpace_discreteRepr('Ephi', typeOfCylinder,rho, q, p, k_0,  m,...
                        a_excitationOfEwave, a_excitationOfHwave));
                    E_z = (k_0/2/pi) * sum(dq_simp.* field__ForBeamDiffraction_InFreeSpace_discreteRepr('Ez', typeOfCylinder,rho, q, p, k_0,  m,...
                        a_excitationOfEwave, a_excitationOfHwave));
                    H_phi = (k_0/2/pi) * sum(dq_simp.* field__ForBeamDiffraction_InFreeSpace_discreteRepr('Hphi', typeOfCylinder,rho, q, p, k_0,  m,...
                        a_excitationOfEwave, a_excitationOfHwave));
                    H_z = (k_0/2/pi) * sum(dq_simp.* field__ForBeamDiffraction_InFreeSpace_discreteRepr('Hz', typeOfCylinder,rho, q, p, k_0,  m,...
                        a_excitationOfEwave, a_excitationOfHwave));F
                    
                S_component = c / (8 * pi) * real(E_phi.* conj(H_z) -  E_z.* conj(H_phi));
                
%             case 'S_phi'
                
            case 'S_z'
                    E_phi = (k_0/2/pi) * sum(dq_simp.* field__ForBeamDiffraction_InFreeSpace_discreteRepr('Ephi', typeOfCylinder, rho, z, q, p, k_0,  m,...
                        a_excitationOfEwave, a_excitationOfHwave));
                    E_rho = (k_0/2/pi) * sum(dq_simp.* field__ForBeamDiffraction_InFreeSpace_discreteRepr('Erho', typeOfCylinder, rho, z, q, p, k_0,  m,...
                        a_excitationOfEwave, a_excitationOfHwave));
                    H_phi = (k_0/2/pi) * sum(dq_simp.* field__ForBeamDiffraction_InFreeSpace_discreteRepr('Hphi', typeOfCylinder, rho, z, q, p, k_0,  m,...
                        a_excitationOfEwave, a_excitationOfHwave));
                    H_rho = (k_0/2/pi) * sum(dq_simp.* field__ForBeamDiffraction_InFreeSpace_discreteRepr('Hrho', typeOfCylinder, rho, z, q, p, k_0,  m,...
                        a_excitationOfEwave, a_excitationOfHwave));
                    
                S_component = c / (8 * pi) * real(E_rho.* conj(H_phi) -  E_phi.* conj(H_rho));
        end
