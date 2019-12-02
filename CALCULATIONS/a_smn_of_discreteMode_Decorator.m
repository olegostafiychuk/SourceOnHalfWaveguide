%%%%% расчёт коэффициента для волн непрерывного спектра в праямом
%%%%% направлении (в сторону среды 2)
%%%%% второе приближение в нахождении коэффициента aplus
function aplus2 = a_smn_of_discreteMode_Decorator(typeOfCylinder, q, p, waveguideParameters, sourceParameters)

m = sourceParameters.m;
d = sourceParameters.d;
j_f = sourceParameters.j_f;
j_z = sourceParameters.j_z;
k_0 = sourceParameters.k_0;
b_0 = sourceParameters.radius;
p_0 = sourceParameters.p_0;

EE1 = waveguideParameters.EEinner;
GG1 = waveguideParameters.GGinner;
HH1 = waveguideParameters.HHinner;
MU1 = waveguideParameters.MUinner;
EE = waveguideParameters.EEouter;
MU = waveguideParameters.MUouter;
a_0 = waveguideParameters.radius;

aplus2 = a_smn_of_discreteMode(typeOfCylinder, q, p_0, p, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, j_f, j_z, d);
