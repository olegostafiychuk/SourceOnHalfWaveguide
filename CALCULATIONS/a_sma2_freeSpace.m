function a_sma = a_sma2_freeSpace(q, p, k_0, p_0, a_0, m, c, I_z, I_f, d)

dJm = (-m)./(k_0* q * a_0).* besselj(-m, k_0* q * a_0) - besselj(-m+1, k_0* q * a_0);
N_p2 = (-1)^(m+2) * c* p./ (k_0.^2 * q);
Ephi2= 1i * dJm;
a_sma = (4 * pi * a_0./ N_p2).* (I_f.* Ephi2).* sin(k_0 * (p - p_0) * d)./ (k_0 * (p - p_0));