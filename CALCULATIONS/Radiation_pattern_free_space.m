clearvars
clc
% вычисление диаграммы направленности несимметричного источника
% typeOfCylinder = 'Gyrotropic';
% typeOfCylinder = 'Isotropic';
systemParameters
% EE1 = 1.0000001;
% GG1 = 0.0000001;
% HH1 = 1.0000001;

step = 1/10000;
teta1 = (0.001:step:0.499999)*pi;
teta2 = (0.500001:step:0.999)*pi;
% teta2 = teta1;
q1 = sin(teta1);
q2 = sin(teta2);
p1 = sqrt(1 - q1.^2);%cos(teta1);
p2 = sqrt(1 - q2.^2);%cos(teta2);%
z  = 0;

a_E_plus  = a_sma1_freeSpace(q1, p1, k_0, p_0, a_0, m, c, j_z, j_f, d);
a_E_minus = a_sma1_freeSpace(q2, -p2, k_0, p_0, a_0, m, c, j_z, j_f, d);
a_H_plus  = a_sma2_freeSpace(q1, p1, k_0, p_0, a_0, m, c, j_z, j_f, d);
a_H_minus = a_sma2_freeSpace(q2, -p2, k_0, p_0, a_0, m, c, j_z, j_f, d);

D_plus  = c/(8*pi*k_0^2).*(cot(teta1)).^2.*(a_E_plus.*conj(a_E_plus) + a_H_plus.*conj(a_H_plus)); % (пока только для
%положительного направления оси z) %
D_minus = c/(8*pi*k_0^2).*(cot(teta2)).^2.*(a_E_minus.*conj(a_E_minus) + a_H_minus.*conj(a_H_minus)); % для отрицательного направления оси z
D_max = max(max(abs(D_plus)),max(abs(D_minus)));

% figure(31)
% %subplot(1,2,1)
% plot(teta1/pi,real(D_plus/D_max)); hold on; %, teta1/pi, imag(D_plus/D_max),'--');  hold on;
% xlabel('teta (\pi rad)'); ylabel('S_r/max(|S_r|)');
% %subplot(1,2,2)
% plot(teta2/pi,real(D_minus/D_max)); hold on; %, teta2/pi, imag(D_minus/D_max),'--'); hold on;
% %xlabel('teta')

teta = [teta1  teta2];
D1 = [D_plus zeros(1,size(teta2,2))];
D2 = [zeros(1,size(teta1,2)) D_minus];
D = (D1 +D2)/D_max;
figure(31)
plot(teta/pi, D,'k'); hold on;
xlabel('teta (\pi rad)'); ylabel('S_r/max(|S_r|)');

% Построение в полярных координатах
teta = [teta1  teta2];
D1 = [D_plus zeros(1,size(teta2,2))];
D2 = [zeros(1,size(teta1,2)) D_minus];
D = (D1 +D2)/D_max;
%D = [D_plus D_minus];
figure(32)
polar(teta, D, 'k'); hold on;
polar(-teta, D, 'k');

% figure
% plot(teta/pi,D)

% figure
% hold on
% plot(teta1/pi,D_plus)
% plot(teta2/pi,D_minus)
% hold off


















