function N_p = Norm_of_continuousWaves(q, p, k_0, c, m, psi, Cm2, Dm2)

%%%%%% norm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_p = (-1)^m * (4 * c / k_0^2) * p./ q.* (Cm2.^2 + Dm2.^2).* psi;
