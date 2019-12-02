function a_sma = a_sma1_freeSpace(q, p, k_0, p_0, a_0, m, c, I_z, I_f, d)

Jm = besselj(-m, k_0* q * a_0);
N_p1 = (-1)^(m+1) * c * p./ (k_0.^2 * q);
Ez  = q.* Jm;
Ephi= -(k_0 * q.^2).^(-1).* ((-p).* (-m) / a_0).* Ez;
a_sma = (4 * pi * a_0./ N_p1).*  (I_z.* Ez + I_f.* Ephi).* sin(k_0 * (p - p_0) * d)./ (k_0 * (p - p_0));