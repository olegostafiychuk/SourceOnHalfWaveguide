function [pp,  Z_val1 ] = GRID_search(p_0, m, EE, GG, HH, k_0, a_0)
%GRID searching
%   поиск корня вручную по сетке
%p_0 = [111 -4]; %[57 -265];  % начальная точка
L = 0.2;%0.005;
tol = 1e-8;%было 1e-10
Z_val = 1;
n = 300;
%while_flag = 0;
n_iter = 0;
while L > tol
        n_iter = n_iter+1;
    X_grid = p_0(1)-L:L/n:p_0(1)+L;
    Y_grid = p_0(2)-L:L/n:p_0(2)+L;
    [X, Y] = meshgrid(X_grid,Y_grid);
    Z = dispeq_gyrotropic_cylinder(X + 1i.*Y, m, EE, GG, HH, k_0, a_0);
    [Z_val1, pos1] = min(min(Z,[],1)); %pos1 - номер столбца, содержащего минимальный элемент
    [Z_val2, pos2] = min(min(Z,[],2)); %pos2 - номер строки, содержащей минимальный элемент
    pp = [X(pos2,pos1) Y(pos2,pos1)];
        if pos1==1 || pos1==size(Z,1) || pos2==1 || pos2==size(Z,2)
            while pos1==1 || pos1==size(Z,1) || pos2==1 || pos2==size(Z,2)%условие нахождения минимума внутри заданной области
                X_grid = pp(1)-L:L/n:pp(1)+L; 
                Y_grid = pp(2)-L:L/n:pp(2)+L;
                [X, Y] = meshgrid(X_grid,Y_grid);
                Z = dispeq_gyrotropic_cylinder(X + 1i.*Y, m, EE, GG, HH, k_0, a_0);
                [Z_val1, pos1] = min(min(Z,[],1)); %pos1 - номер столбца, содержащего минимальный элемент
                [Z_val2, pos2] = min(min(Z,[],2)); %pos2 - номер строки, содержащей минимальный элемент
                pp = [X(pos2,pos1) Y(pos2,pos1)];
            end
        end
    
    p_0 = pp;
    L = L/(0.5*n);
        if Z_val1 == Z_val || L == 0 || isnan(Z_val)
            %while_flag = 1;
            %disp('Невозможно найти // заданный минимум функции// при заданной точности');
            break
        end
%     Z_val = Z_val1;
    Z_val = (Z_val1 + Z_val2)./2;
end
% p_err = 0;
% switch while_flag
%     case 1
%         p_err = 1;
%         disp('Невозможно найти ноль функции. Выберите другую сетку');
%         disp('Найденное минимальное значение функции внутри заданной области:');
%         result_fmin = ['f_min = ',num2str(Z_val)];
%         disp(result_fmin);
%         disp('Количество итераций:');
%         NN = ['n = ',num2str(n_iter)];
%         disp(NN);
%     otherwise
%         result_w = ['w=',num2str(w0/wLH),'wLH'];
%         result_p = ['p = ',num2str(p)];
%         result_fmin = ['f_min = ',num2str(Z_val)];
%         PP = ['p""/p" = ',num2str(imag(p)/real(p))];
%         disp(result_w);
%         disp(result_p);
%         disp(result_fmin);
%         disp(PP);
%         disp('Количество итераций:');
%         NN = ['n = ',num2str(n_iter)];
%         disp(NN);
% %     format long
% end



end

