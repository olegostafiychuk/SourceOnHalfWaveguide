% /*===========================================================================
% 
% DESCRIPTION
%       ѕрограмма вычисл€ет зависимость "диспесионной" функции от
%       действительной и мнимой частей посто€нной распространени€ при
%       фиксированной частоте.
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

clearvars
clc
% global EE GG HH k_0 a_0 m

m = 1;


eps = 1;
y0  = 0;
x0  = 1;
xmin =  x0 - eps;
xmax =  x0 + eps;
ymin =  y0 - eps;
ymax =  y0 + eps;
Npntx = 100;
Npnty = 100;


systemParameters

EE = EE1
GG = GG1
HH = HH1
R   = k_0 * a_0;

% %%refractive-index surface%%%%%%%%%%%%%%%%%%%%%%
% p = 1:0.1:50;
% mainq = EE.^2 - GG.^2 + EE.* HH - (HH + EE).* p.^2;
%     radq = sqrt((HH - EE).^2 * p.^4 + 2 * ((GG.^2).* (HH + EE) - EE.* (HH - EE).^2) * p.^2 +...
%         (EE.^2 - GG.^2 - EE.* HH).^2);
%     q1 = sqrt(0.5 * (mainq - radq)./ EE);
%     q2 = sqrt(0.5 * (mainq + radq)./ EE);
% plot(p,q1,'r',p,q2,'b');
  
    

[ff1, fval] = fminsearch(@(x) dispeq_gyrotropic_cylinder(x(1) + 1i * x(2), m, EE, GG, HH, k_0, a_0), [0.366 0], optimset('TolX',1e-8))
%fval

[px,py] = meshgrid(xmin : (xmax - xmin) / Npntx :xmax, ymin : (ymax - ymin) / Npnty :ymax);
%         zz1 = (dispeqh1_new(px + 1i * py));
p = px + 1i * py;
%
% [qx,qy] = meshgrid(xmin : (xmax - xmin) / Npntx :xmax, ymin : (ymax - ymin) / Npnty :ymax);
% p = sqrt(1-(qx+1i*qy).^2);

%         zz1 = (dispeq_isotropic_cylinder(px + 1i * py));
        zz1 = (dispeq_gyrotropic_cylinder(p, m, EE, GG, HH, k_0, a_0));
%         zz1 = ((GG1^2 - (p.^2 - EE1).^2));
        zz2 = abs(zz1);
        
        
%         zz10 = (dispeqh1_TE_new(px + i * py, 1, 0));
%         zz01 = (dispeqh1_TE_new(px + i * py, 0, 1));
%         zz00 = (dispeqh1_TE_new(px + i * py, 0, 0)); 

% for I = 1 : Npntx
%     for J = 1 : Npnty
% %         zz11(J,I) = (dispeqh1_TE([px(J, I) py(J, I)], 1, 1));
%         zz10(J,I) = (dispeqh1_TE([px(J, I) py(J, I)], 1, 0));
% %         zz01(J,I) = (dispeqh1_TE([px(J, I) py(J, I)], 0, 1));
% %         zz00(J,I) = (dispeqh1_TE([px(J, I) py(J, I)], 0, 0));
%     end
% end

% zz2 = log(zz2);
figure
mesh(px,py, zz2)
title('zz11');
xlabel('Re p');
ylabel('Im p');
set(gca,'ZScale','log');
view(30,30)

%%
%%%%% search of mode constant
%%% the search of transverse constant of surface modes
p_min = 1;
p_Max = 2;
    p = [p_min:0.1:p_Max];
    p_n = zeros(size(p,2),2);
    p_n_pred = zeros(size(p,2),1);
    fval = zeros(size(p,2),1);
    delta = 0.1;%1.5;
for in = 1:size(p,2)
%     [p_n(in,:), fval(in)] = fminsearch(@(x) dispeq_gyrotropic_cylinder(x(1) + 1i * x(2), m, EE, GG, HH, k_0, a_0),...
%                 [p(in),0], optimset('TolX',1e-8));
    [p_n(in,:), fval(in)] = fminbnd(@(x) dispeq_gyrotropic_cylinder(x, m, EE, GG, HH, k_0, a_0),...
                p(in) - delta, p(in) + delta, optimset('TolX',1e-8));
            
%     [p_n(in,:), fval(in)] = fminsearch(@(x) dispeq_gyrotropic_cylinder(x(1) + 1i * x(2), m, EE, GG, HH, k_0, a_0),...
%                p_n(in,:), optimset('TolX',1e-8));
end
% 
%%% sorting of transverse propagation constant
sizep_n=size(p_n);
for i=1:sizep_n(1)-1
for j=i+1:sizep_n(1)
    if p_n(i,1)>p_n(j,1)
        p_ntemp=p_n(i,:);
        p_n(i,:)=p_n(j,:);
        p_n(j,:)=p_ntemp;
    end
end
end
% 
% %%% forming no equal mode 
zz = p_n(1,:);
sizey = size(p_n(:,1));
for t = 1:sizey(1)-1
    sizezz = size(zz(:,1));
          y1 = p_n(t,1);
          y2 = p_n(t + 1,1);
          if abs(y1-y2)>0.0001 
             zz(sizezz(1) + 1,1) = p_n(t+1,1);
%              zz(sizezz + 1,2) = q_n(t,2);
          end
end
p_n = zz;
%p_n = p_n(:,p_n(:,1)<p_Max);
p_n = p_n(p_n(:,1)<p_Max,:);
% 
% 
% 
% % %%% the search of transverse constant of complex modes
% % p_min = 0;
% % p_Max = 15;
% %     p = [p_min:0.1:p_Max];
% %     
% % for in = 1:size(p,2)
% %     p_n(in,:) = fminsearch(@(x) dispeq_gyrotropic_cylinder(x(1) + 1i * x(2)),...
% %                 [1, p(in)], optimset('TolX',1e-8))
% % end
% % 
% % %%% sorting of transverse propagation constant
% % sizep_n=size(p_n);
% % for i=1:sizep_n(1)-1
% % for j=i+1:sizep_n(1)
% %     if p_n(i,1)>p_n(j,1)
% %         p_ntemp=p_n(i,:);
% %         p_n(i,:)=p_n(j,:);
% %         p_n(j,:)=p_ntemp;
% %     end
% % end
% % end
% % 
% % %%% forming no equal mode 
% % zz = p_n(1,:);
% % sizey = size(p_n(:,1));
% % for t = 1:sizey(1)-1
% %     sizezz = size(zz(:,1));
% %           y1 = p_n(t,1);
% %           y2 = p_n(t + 1,1);
% %           if abs(y1-y2)>0.00001 
% %              zz(sizezz(1) + 1,:) = p_n(t+1,:);
% % %              zz(sizezz + 1,2) = q_n(t,2);
% %           end
% % end
% % p_n = zz



