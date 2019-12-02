function component = field_of_discreteMode_Decorator(componentOfField, typeOfCylinder, r, z, q, p, a_p_field, waveguideParameters, sourceParameters)
EE1 = waveguideParameters.EEinner;
GG1 = waveguideParameters.GGinner;
HH1 = waveguideParameters.HHinner;
MU1 = waveguideParameters.MUinner;
EE = waveguideParameters.EEouter;
MU = waveguideParameters.MUouter;
a_0 = waveguideParameters.radius;

m = sourceParameters.m;
k_0 = sourceParameters.k_0;

component = field_of_discreteMode(componentOfField, typeOfCylinder, r, q, p, a_p_field, k_0, a_0, EE1, GG1, HH1, MU1, EE, MU, m, z);
  

     
     





