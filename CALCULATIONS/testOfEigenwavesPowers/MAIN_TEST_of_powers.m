clear all
tic
systemParameters

p_n__of_descreteMode_of_gyrotropicCyl;

P_n = zeros(size(p_n,1),1);
innerMethod = zeros(size(p_n,1),1);
outerMethod = zeros(size(p_n,1),1);

P_n_withLoss = zeros(size(p_n,1),1);


for in = 1:size(p_n,1)
    [Pn, innerAnalytical, outerAnalytical] = powerOfDiscreteSpectrumMode(p_n(in), waveguideParameters, sourceParameters);
    P_n(in) = Pn;
    innerMethod(in) = innerAnalytical;
    outerMethod(in) = outerAnalytical;

    [Pn, innerAnalytical, outerAnalytical] = powerOfDiscreteSpectrumMode(p_n(in)* (1 - 1i*1e-13), waveguideParameters, sourceParameters);
    P_n_withLoss(in) = Pn;
end

figure(1)
hold on
plot(p_n, P_n);
plot(p_n, P_n_withLoss, 'ro');
hold off

figure(2)
hold on
plot(p_n, innerMethod, 'ro');
plot(p_n, outerMethod, 'bo');
hold off

toc
