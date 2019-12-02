%%%%% расчёт коэффициента для волн непрерывного спектра в прямом
%%%%% направлении (в сторону среды 2)
%%%%% второе приближение в нахождении коэффициента aplus
function component = field__ContinuesWaves_discreteRepr(componentOfField, typeOfCylinder,r, z, q, p,...
    a_p_field_1_forw, a_p_field_1_back, a_p_field_2_forw, a_p_field_2_back,...
    waveguideParameters, sourceParameters, coeffsOfField)

m = sourceParameters.m;
k_0 = sourceParameters.k_0;
k = k_0;

EE1 = waveguideParameters.EEinner;
GG1 = waveguideParameters.GGinner;
HH1 = waveguideParameters.HHinner;
MU1 = waveguideParameters.MUinner;
EE = waveguideParameters.EEouter;
MU = waveguideParameters.MUouter;
a_0 = waveguideParameters.radius;

psi_forward_1 = coeffsOfField.psi_forward_1;
psi_backward_1 = coeffsOfField.psi_backward_1;
psi_forward_1_transp = coeffsOfField.psi_forward_1_transp;
psi_backward_1_transp = coeffsOfField.psi_backward_1_transp;
psi_forward_2 = coeffsOfField.psi_forward_2;
psi_backward_2 = coeffsOfField.psi_backward_2;
psi_forward_2_transp = coeffsOfField.psi_forward_2_transp;
psi_backward_2_transp = coeffsOfField.psi_backward_2_transp;
B_1_forward_1 = coeffsOfField.B_1_forward_1;
B_2_forward_1 = coeffsOfField.B_2_forward_1;
Cm2_forward_1 = coeffsOfField.Cm2_forward_1;
Dm2_forward_1 = coeffsOfField.Dm2_forward_1;
B_1_backward_1 = coeffsOfField.B_1_backward_1;
B_2_backward_1 = coeffsOfField.B_2_backward_1;
Cm2_backward_1 = coeffsOfField.Cm2_backward_1;
Dm2_backward_1 = coeffsOfField.Dm2_backward_1;
B_1_forward_2 = coeffsOfField.B_1_forward_2;
B_2_forward_2 = coeffsOfField.B_2_forward_2;
Cm2_forward_2 = coeffsOfField.Cm2_forward_2;
Dm2_forward_2 = coeffsOfField.Dm2_forward_2;
B_1_backward_2 = coeffsOfField.B_1_backward_2;
B_2_backward_2 = coeffsOfField.B_2_backward_2;
Cm2_backward_2 = coeffsOfField.Cm2_backward_2;
Dm2_backward_2 = coeffsOfField.Dm2_backward_2;

component = field__ForBesselBeamRepresentation_Continues_discreteRepr(componentOfField, typeOfCylinder,r, q, 0, k_0, k, a_0, EE1, GG1, HH1, MU1, 0, 0, EE, MU, m, z, 0, 0,...
          a_p_field_1_forw, a_p_field_1_back, a_p_field_2_forw, a_p_field_2_back,...
          psi_forward_1, psi_backward_1, psi_forward_1_transp, psi_backward_1_transp,...
          psi_forward_2, psi_backward_2, psi_forward_2_transp, psi_backward_2_transp,...
          B_1_forward_1,B_2_forward_1,Cm2_forward_1,  Dm2_forward_1,...
          B_1_backward_1,B_2_backward_1,Cm2_backward_1, Dm2_backward_1,...
          B_1_forward_2,B_2_forward_2,Cm2_forward_2,  Dm2_forward_2,...
          B_1_backward_2,B_2_backward_2,Cm2_backward_2, Dm2_backward_2);
