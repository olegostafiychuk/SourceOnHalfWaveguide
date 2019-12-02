%%source resistance in free space and in the presence of the waveguide on frequency

function [R_alp1, R_alp2, R_E_freespace, R_H_freespace] = sourceResistance_in_presence_waveguide_and_freeSpace(w0s, typeOfCylinder,...
    H0, w_H, typeOfmedia, dd, EE_0, p_0, a_0, MU1, EE, MU, c, m, j_f, j_z, I_0, d, toOhms)

    
I = 0;
for D_Q = w0s

w_0 = D_Q;
k_0 = w_0 / c;
R = a_0 * k_0;
        
% %%%%% EE при учете потерь
[EE1, GG1, HH1, c] = channelparameters_sources(H0, typeOfmedia, w_0, dd, EE_0);

    I=I+1;
    Qpts(I) = w_0 / w_H;
    
    %%% ищем моды с помощью матлабовской функции fminserch
    
    delta = 0;
    pp = @(q) sqrt(1 - q.^2);
    N_p1 = @(q, p) (-1)^(m+1) * c * p./ (k_0.^2 * q);
    N_p2 = @(q, p) (-1)^(m+2) * c * p./ (k_0.^2 * q);

    %%%%% for free space
    P_Ewave_plus  =  -1/4 * quadgk(@(q) abs(a_sma1_freeSpace(q, pp(q), k_0, p_0, a_0, m, c, j_z, j_f, d)).^2.*...
        N_p1(q, pp(q)), delta, 1);
    P_Ewave_minus  =   1/4 * quadgk(@(q) abs(a_sma1_freeSpace(q,-pp(q), k_0, p_0, a_0, m, c, j_z, j_f, d)).^2.*...
        N_p1(q,-pp(q)), delta, 1);

    P_Hwave_plus  =  1/4 * quadgk(@(q) abs(a_sma2_freeSpace(q, pp(q), k_0, p_0, a_0, m, c, j_z, j_f, d)).^2.*...
        N_p2(q, pp(q)), delta, 1);
    P_Hwave_minus  =  -1/4 * quadgk(@(q) abs(a_sma2_freeSpace(q,-pp(q), k_0, p_0, a_0, m, c, j_z, j_f, d)).^2.*...
        N_p2(q,-pp(q)), delta, 1);

    %%%%% for waveguide
    P_alp1_forw  =  quadgk(@(q) P_of_continuousWaves_alpha1(typeOfCylinder,...
        q, p_0, pp(q), k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, 0, j_f, j_z, d), delta, 1);
    P_alp1_back  = -quadgk(@(q) P_of_continuousWaves_alpha1(typeOfCylinder,...
        q, p_0,-pp(q), k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, 0, j_f, j_z, d), delta, 1);
    
    P_alp2_forw  = -quadgk(@(q) P_of_continuousWaves_alpha2(typeOfCylinder,...
        q, p_0, pp(q), k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, 0, j_f, j_z, d), delta, 1);
    P_alp2_back  = quadgk(@(q) P_of_continuousWaves_alpha2(typeOfCylinder,...
        q, p_0,-pp(q), k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, c, m, 0, j_f, j_z, d), delta, 1);
    
    R_E_freespace(I) = abs((P_Ewave_plus + P_Ewave_minus)) * 2 / I_0^2 * toOhms;
    R_alp1(I)        = abs((P_alp1_forw + P_alp1_back))  * 2 / I_0^2   * toOhms;
    R_H_freespace(I) = abs((P_Hwave_plus + P_Hwave_minus)) * 2 / I_0^2 * toOhms;
    R_alp2(I)        = abs((P_alp2_forw + P_alp2_back))  * 2 / I_0^2   * toOhms;
end
