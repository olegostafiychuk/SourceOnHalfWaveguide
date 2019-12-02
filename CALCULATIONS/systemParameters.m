typeOfmedia = 'Plasma';
% typeOfmedia = 'Graphen';
w = 0;
channelparameters;

%%% parameters of the source
w_0 = 0.025 * w_H;
% w_0 = 0.022 * w_H;
% w_0    = 0.01 * w_H;
% w_0    = 0.002 * w_H;
% p_0 = 0.5;
% p_0 = 31.3317352095852;
% p_0 = 27.4089166804486;
p_0 = 0;
% p_0 = 200.7;
m = 1;
k_0 = w_0 / c;
% d = 4 * a_0;
d = 2 * pi / k_0 / 100;
j_f = 1e7 / (2 * d) * 0;
j_z = 2e6 / (2 * pi * a_0);
I_0 = sqrt((j_z * 2 * pi * a_0)^2 + (j_f * 2 * d)^2);
L     = - 2 * pi / k_0 * 2; %%% z coordinate of center of source

sourceParameters.m = m;
sourceParameters.d = d;
sourceParameters.j_f = j_f;
sourceParameters.j_z = j_z;
sourceParameters.p_0 = p_0;
sourceParameters.I_0 = I_0;
sourceParameters.w_0 = w_0;
sourceParameters.k_0 = k_0;
sourceParameters.radius = a_0;
sourceParameters.zCoordinate = L;

%%% paramenter of the outer media
EE = 1;
MU = 1;

%%% parameters of the waveguide media
H0 = 800;
[EE1, GG1, HH1, c] = channelparameters_sources(H0, typeOfmedia, w_0, 0, 0);
MU1 = 1;
% a_0 = 0.001;

waveguideParameters.EEinner = EE1;
waveguideParameters.GGinner = GG1;
waveguideParameters.HHinner = HH1;
waveguideParameters.MUinner = MU1;
waveguideParameters.typeOfmedia = typeOfmedia;
waveguideParameters.EEouter = EE;
waveguideParameters.MUouter = MU;
waveguideParameters.radius = a_0;

%%% paramenter of right media
EE2 = 1;
MU2 = 1;
