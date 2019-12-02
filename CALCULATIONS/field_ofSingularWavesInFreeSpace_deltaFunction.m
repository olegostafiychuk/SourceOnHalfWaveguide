function component = field_ofSingularWavesInFreeSpace_deltaFunction(componentOfField, r, q, p, k_0, m, z, a_p_field_Ewave_forw, a_p_field_Hwave_forw)

        %%%%%%%%%%%%%%%%%%%%%%% E-polarized beam %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Ez_inc    = a_p_field_Ewave_forw.* q.*(besselj(m, k_0.* r * q)); 
  dEz_q_inc = a_p_field_Ewave_forw.* k_0.* q.*q.*((besselj(m,   k_0.* r* q) * (m))./ (k_0.* r* q)  - besselj(m + 1,   k_0.* r* q));

%%%%%%%%%%%%%%%%%%%%%%% H-polarized beam %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Hz_inc    = a_p_field_Hwave_forw.* q.*(besselj(m, k_0.* r * q));
  dHz_q_inc = a_p_field_Hwave_forw.* k_0.* q.* q.*((besselj(m,   k_0.* r* q) * (m))./ (k_0.* r* q)  - besselj(m + 1,   k_0.* r* q));
  
  A0 = 1./ (k_0 * (1 - p.^2));
        switch(componentOfField)
            case 'Ez'
                component = Ez_inc;
            case 'Hz'
                component = Hz_inc;
            case 'Ephi'
                component = A0.* ((-p.* (m./r)).* Ez_inc + 1i * dHz_q_inc);
            case 'Hphi'
                component = A0.* (-1i * dEz_q_inc - (p.* (m./r)).* Hz_inc);
            case 'Erho'
                component = A0.* (-1i * p.* dEz_q_inc - (m./r).* Hz_inc);
            case 'Hrho'
                component = A0.* ((m./r).* Ez_inc - 1i * p.* dHz_q_inc);
        end

component = k_0/2/pi * component.* exp(-1i * k_0 * p * z);







