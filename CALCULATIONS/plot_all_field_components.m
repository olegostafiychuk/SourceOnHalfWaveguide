function ff = plot_all_field_components(figNumber, rho, a_0, sgn, wid, Ez_inc, Ephi_inc, Erho_inc,...
                                        Hz_inc, Hphi_inc, Hrho_inc)

% wid = 1;                                    
%%%%%%%%%%%%%%%%%%%%% painting of field %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(figNumber)
% subplot(2,3,1);
figure(231);
hold on
% axes('FontSize',16)
plot(rho/a_0, real(Ez_inc),['b' sgn],'LineWidth',wid)
plot(rho/a_0, imag(Ez_inc),['r' sgn] ,'LineWidth',wid)
hold off
title('Ez')
box on
% subplot(2,3,2);
figure(232);
hold on
% axes('FontSize',16)
plot(rho/a_0, real(Ephi_inc),['b' sgn],'LineWidth',wid)
plot(rho/a_0, imag(Ephi_inc),['r' sgn] ,'LineWidth',wid)
hold off
title('Ephi')
box on
% subplot(2,3,3);
figure(233);
hold on
% axes('FontSize',16)
plot(rho/a_0, real(Erho_inc),['b' sgn],'LineWidth',wid)
plot(rho/a_0, imag(Erho_inc),['r' sgn] ,'LineWidth',wid)
hold off
title('Erho')
box on
% subplot(2,3,4);
figure(234);
hold on
% axes('FontSize',16)
plot(rho/a_0, real(Hz_inc),['b' sgn],'LineWidth',wid)
plot(rho/a_0, imag(Hz_inc),['r' sgn] ,'LineWidth',wid)
hold off
title('Hz')
box on
% subplot(2,3,5);
figure(235);
hold on
% axes('FontSize',16)
plot(rho/a_0, real(Hphi_inc),['b' sgn],'LineWidth',wid)
plot(rho/a_0, imag(Hphi_inc),['r' sgn] ,'LineWidth',wid)
hold off
title('Hphi')
box on
% subplot(2,3,6);
figure(236);
hold on
% axes('FontSize',16)
plot(rho/a_0, real(Hrho_inc),['b' sgn],'LineWidth',wid)
plot(rho/a_0, imag(Hrho_inc),['r' sgn] ,'LineWidth',wid)
hold off
title('Hrho')
box on
