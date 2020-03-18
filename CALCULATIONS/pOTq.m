% duct parameters:
global m R EE GG HH ee Q n Q1 Q2 n1 n2 alp alp1 alp2
% m=1; D_Q=0.101; D_k=5.623; D_d=1.2; D_g=(4.3e-3)^2; 
% D_R=2.344; V_wH=0.00; v_wH=0.00;
%

systemParameters
m = 1
EE = EE1
GG = GG1
HH = HH1


q=[0:0.001:50];

Rq = sqrt(0.25*(1 - EE./HH)^2.*q.^4 - GG^2/HH.*q.^2 + GG^2);
po = sqrt(EE - q.^2.* 0.5.* (1 + EE./HH) - Rq);
px = sqrt(EE - q.^2.* 0.5.* (1 + EE./HH) + Rq);

figure(11)
hold on
plot(q, real(px), 'blue', q, imag(px), 'r')
hold off

figure(12)
hold on
plot(q, real(po), 'blue', q, imag(po), 'r')
hold off

figure(17)
hold on
plot(q, real(Rq), 'blue', q, imag(Rq), 'r')
hold off





p = [-10:0.001:10];

mainq = EE.^2 - GG.^2 + EE.* HH - (HH + EE).* p.^2; 
radq = sqrt((HH - EE).^2 * p.^4 + 2 * ((GG.^2).* (HH + EE) - EE.* (HH - EE).^2) * p.^2 +...
    (EE.^2 - GG.^2 - EE.* HH).^2);
    
q1 = sqrt(0.5 * (mainq - radq)./ EE);
q2 = sqrt(0.5 * (mainq + radq)./ EE);
% q1 = q1.* (2*(imag(q1) <= 0)-1);
% q2 = q2.* (2*(imag(q2) <= 0)-1);

n1 = -(EE.* (p.* GG).^(-1)).* (p.^2 + q1.^2 + (GG.^2)./ EE - EE);
n2 = -(EE.* (p * GG).^(-1)).* (p.^2 + q2.^2 + (GG.^2)./ EE - EE);

figure(13)
hold on
plot(real(q1),p, 'blue', imag(q1),p, 'r')
hold off

figure(14)
hold on
plot(real(q2), p, 'blue', imag(q2), p,'r')
hold off



