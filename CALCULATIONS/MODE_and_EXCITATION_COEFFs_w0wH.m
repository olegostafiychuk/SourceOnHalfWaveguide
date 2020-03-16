% /*===========================================================================
% 
% DESCRIPTION
%       Программа вычисляет зависимости постоянных распространения от частоты 
%       и строит их графики. Для затравки задаются значения постоянных
%       распространения при конкретной частоте.
% 
%   Copyright (c) 2005-2013 by Vasiliy Es'kin. All Rights Reserved.
% ===========================================================================*/
% 
%                       EDIT HISTORY FOR FILE
% 
%   This section contains comments describing changes made to the module.
%   Notice that changes are listed in reverse chronological order.
% 
% when       who              what, where, why
% --------   ---       ----------------------------------------------------------
% 11/09/2007 Vasiliy Es'kin   Create programma.
% ==========================================================================*/


%вычисление дисперсионных кривых
clearvars     %%% очищаем память перед новыми вычислениями
clc

tic

% global m R EE GG HH a_0 w_H w_p O_H O_p c k_0
systemParameters

%R   = k_0 * a_0;


% upper_Bound = 4000;
% N = 400;%
% [dq_simp, q_back, p] = quadratureMethod_forIntegralEqs('simpson', N, -1i, 1e-8, upper_Bound);
% p_back = sqrt(1 - q_back.^2);
% p_back = real(p_back) - 1i * abs(imag(p_back));
%%%%% the coefficient of forward waves in waveguend space
% q = q_back;
% p = p_back;


typeOfCylinder = 'Gyrotropic';
%%%%% number of scattering mode
% Nn = 7;
% Nn = 35;

p_n__of_descreteMode_of_gyrotropicCyl_NEW; % бесстолкновительный случай
% p_n = p_n(real(p_n)< 25.5);
% p_n = p_n((150<real(p_n))&(real(p_n)< 160));
%%% задаём постоянные распространения мод, найденные для
%%% бесстолкновительной плазмы

% p_n_losses;

x0(:,1) = p_n; %p_n_real(:, end);
x0(:,2) = 0; %p_n_imag(:, end);%- 1e-15;

% pp(:,1) = p_n;
% pp(:,2) = 0;%- 1e-15;

%w_Gran_wH = sqrt((w_p^2 + w_H^2)/2)/w_H;

clear Qpts p_ p__ Po Pi

%%% задаём начальное и конечное значение частоты
Qmin = 0;%Nu_e_vec(end);
Qmax = 1.5e-2 * w_H;%5e-3 * w_H;%
n = 150;
% Qpts=[Qmin:0.01*(Qmax-Qmin):Qmax];
deltax = 0.1;
I=0;
for D_Q = [Qmin:(Qmax-Qmin)/n:Qmax];

Nu_e = D_Q;
k_0 = w_0 / c;
%R = a_0 * k_0;
        
% %%%%% EE при учете потерь
Nu_i = 0;
EE = 1 + w_p.^2 * (w_0 - 1i*Nu_e)./ ((w_H.^2 - (w_0 - 1i*Nu_e).^2).* w_0) +...
         O_p.^2 * (w_0 - 1i*Nu_i)./ ((O_H.^2 - (w_0 - 1i*Nu_i).^2).* w_0);
        
GG =  - w_p.^2 * w_H./ ((w_H.^2 - (w_0 - 1i*Nu_e).^2).* w_0) +...
              O_p.^2 * O_H./ ((O_H.^2 - (w_0).^2).* w_0);
          
HH = 1 - w_p.^2./ ((w_0 - 1i*Nu_e).* w_0) - O_p.^2./ ((w_0).* w_0);

%[EE, GG, HH, c] = channelparameters_sources(H0, typeOfmedia, w_0, 0, 0);

    I=I+1
    Qpts(I) = Nu_e;
    Pi(I) = (sqrt(EE-GG));
    Pc(I) = sqrt(EE - ((EE + HH).* GG.^2).* (EE - HH).^(-2)  +...
         2 * ((EE - HH).^(-2)).* sqrt((EE.* HH.* GG.^2).* (GG.^2 - (EE - HH).^2)));
    
     
           dNu_e = (Qmax - Qmin)/n;
     
    %%% ищем моды с помощью матлабовской функции
    for J = 1:size(x0,1)
%         x0(J,:)=fminsearch(@(x) (dispeq_gyrotropic_cylinder(x)),[x0(J,:)],optimset('TolX',1e-12));
          %x0(J,:) = fminbnd(@(x) (dispeq_gyrotropic_cylinder(x)),x0(J,:)-deltax,x0(J,:)+deltax, optimset('TolX',1e-8));
        
%         x0(J,:) = fminbnd(@(x) dispeq_gyrotropic_cylinder(x, m, EE, GG, HH, k_0, a_0),...
%                     x0(J,:)-deltax,x0(J,:)+deltax, optimset('TolX',1e-8));
%         x0(J,:) = fminsearch(@(x) dispeq_gyrotropic_cylinder(x(1)+1i*x(2), m, EE, GG, HH, k_0, a_0), x0(J,:), optimset('TolX',1e-8));
        
        %%%%handle grid search%%%%%%
%         [pp(J,:),  Z_val1 ] = GRID_search(pp(J,:), m, EE, GG, HH, k_0, a_0);
        %pp0(J,:) = pp(J,:);
%         x0(J,:) = fminsearch(@(x) dispeq_gyrotropic_cylinder(x(1)+1i*x(2), m, EE, GG, HH, k_0, a_0), pp(J,:), optimset('TolX',1e-8));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%fminsearchbnd  -5*Nu_e/w_0
        LB = [x0(J,1)-0.1 x0(J,2)-0.1-p_n(J,1)*dNu_e/w_0];%LB = [x0(J,1)-0.1-5*Nu_e/w_0 x0(J,2)-0.1-x0(J,1)*Nu_e/w_0];
        UB = [x0(J,1)+0.1 x0(J,2)+0.1+p_n(J,1)*dNu_e/w_0];%UB = [x0(J,1)+0.1+5*Nu_e/w_0 x0(J,2)+0.1+x0(J,1)*Nu_e/w_0];
        x0(J,:) = fminsearchbnd(@(x) dispeq_gyrotropic_cylinder(x(1)+1i*x(2), m, EE, GG, HH, k_0, a_0),...
            x0(J,:), LB, UB, optimset('TolX',1e-6));        
        
        %pn_real(J,I) = pp(J,1);
        %pn_imag(J,I) = pp(J,2);
        
        p_(J,I)=x0(J,1);
        p__(J,I)=x0(J,2);
    end
    
          
end
toc

figure(2)
hold on
for J=1:size(p_,1)
      plot(Qpts/w_H,p_(J,:),'-*'); hold on
      text(Qpts(1), p_(J,1),['n=' num2str(J)]);
      plot(Qpts/w_H,real(Pi),'k--')
end

figure(3)
hold on
for J=1:size(p_,1)
      plot(Qpts/w_H,p__(J,:),'-*')
      text(Qpts(151)/w_H, p__(J,end),['n=' num2str(J)]);
      plot(Qpts/w_H,imag(Pi),'k--')
end

figure(4)
plot(p_(:,end),p__(:,end),'*'); hold on
plot(0:1:1000,-[0:1:1000]*Qpts(end)/w_0,'k--');
plot(0:1:1000,-[0:1:1000]*Qpts(end)/w_H,'k--');
% hold off

% figure(4)
% hold on
% for J=1:size(p_,1)
%    NofmaxB = 1;
% %     NofmaxB = (abs(real(Bn_back_n(J,:))) == max(abs(real(Bn_back_n(J,:)))));
% %     plot(Qpts,p_(J,:),'b-')
%       plot(Qpts,abs(Bn_back_n(J,:)))
%       text(Qpts(NofmaxB), abs(Bn_back_n(J,NofmaxB))* 1.03,['n=' num2str(J)]);
% end
% hold off


% subplot(2,1,2)
% figure
% % plot(Qpts,imag(Pi),'g--')
% hold on
% for J=1:size(p__,1)
% %     plot(Qpts,-p__(J,:))
%     plot(Qpts,-p__(J,:),'black-')
% end
% hold off


% set(gca,'YScale','log');
% subplot(3,1,3)
% plot(Qpts,real(Q1),'r-',Qpts,imag(Q1),'g')

toc
