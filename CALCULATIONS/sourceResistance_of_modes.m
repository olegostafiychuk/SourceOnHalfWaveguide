%%source resistance in free space and in the presence of the waveguide on frequency

function [R_of_discreteMode, p_n_new, Pin, Pcin] = sourceResistance_of_modes(w0s, p_n, typeOfCylinder,...
    H0, w_H, typeOfmedia, dd, EE_0, p_0, a_0, c, m, j_f, j_z, I_0, d, toOhms)

% global m R EE GG HH k_0 a_0

R_of_discreteMode = zeros(size(p_n,1),size(w0s,2));

deltax = 0.1;
I = 0;
for D_Q = w0s


w_0 = D_Q;
k_0 = w_0 / c;
R = a_0 * k_0;
x0 = p_n;
        
% %%%%% EE при учете потерь
[EE1, GG1, HH1, c] = channelparameters_sources(H0, typeOfmedia, w_0, dd, EE_0);
EE = EE1;
GG = GG1;
HH = HH1;

    I=I+1;
%     Qpts(I) = w_0 / w_H;
    
    Pin(I)  = (sqrt(EE-GG));
    Pcin(I) = sqrt(EE - ((EE + HH).* GG.^2).* (EE - HH).^(-2)  +...
         2 * ((EE - HH).^(-2)).* sqrt((EE.* HH.* GG.^2).* (GG.^2 - (EE - HH).^2)));
    
    %%% ищем моды с помощью матлабовской функции fminserch
    
    %%% ищем моды с помощью матлабовской функции fminserch
    for J=1:size(x0,1)
        deltax = 0.1;
%         x0(J,:)=fminsearch(@(x) (dispeq_gyrotropic_cylinder(x)),[x0(J,:)],optimset('TolX',1e-12));
        x0(J,:)=fminbnd(@(x) (dispeq_gyrotropic_cylinder(x, m, EE, GG, HH, k_0, a_0)),x0(J,:)-deltax,x0(J,:)+deltax, optimset('TolX',1e-8));
        p_n_new(J,I)=x0(J,1);
        
%         while((I>1 && abs(p_n_new(J,I-1)- Pin(I-1))<deltax && abs(p_n_new(J,I)- Pin(I))<deltax) || (J>1 && abs(p_n_new(J,I)- p_n_new(J-1,I))<deltax))
%             deltax = deltax/10;
%             x0(J,:)=fminbnd(@(x) (dispeq_gyrotropic_cylinder(x, m, EE, GG, HH, k_0, a_0)),p_n_new(J,I-1)-deltax,p_n_new(J,I-1)+deltax, optimset('TolX',1e-8));
%             p_n_new(J,I)=x0(J,1);
%         end
    end
    
    p_n = x0(:);
    q_n = sqrt(1-p_n.^2);
    q_n = q_n.* (2*(imag(q_n) <= 0)-1);
    
    for in = 1:size(q_n,1)
        a_smn_forw = a_smn_of_discreteMode(typeOfCylinder, q_n(in), p_0,  p_n(in), k_0, a_0, EE1, GG1, HH1, 1, 1, 1, m, j_f, j_z, d);
        a_smn_back = a_smn_of_discreteMode(typeOfCylinder, q_n(in), p_0, -p_n(in), k_0, a_0, EE1, GG1, HH1, 1, 1, 1, m, j_f, j_z, d);
        
        R_of_discreteMode(in,I) =(P_of_discreteMode(typeOfCylinder, q_n(in), p_n(in), k_0, a_0, EE1, GG1, HH1, 1, 1, 1, m, a_smn_forw)+...
                                 -P_of_discreteMode(typeOfCylinder, q_n(in),-p_n(in), k_0, a_0, EE1, GG1, HH1, 1, 1, 1, m, a_smn_back)).*...
                                  2 / I_0^2 * toOhms;
    end
end
