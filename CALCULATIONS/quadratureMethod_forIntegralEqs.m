function [dq_simp, q_back, p] = quadratureMethod_forIntegralEqs(method, N, q_0, q_start, upper_Bound)


switch(method)
    case 'trapeciech'
    %%%%%%%%%% the formula of trapiech    
    if(q_0/1i ~= imag(q_0) && q_0 > q_start)
        dq = (q_0-q_start)/N;
        q_back1 = [q_start:dq:q_0-dq];
        dq2 = (2-q_0-dq)/N;
        q_back2 = [q_0+dq:dq2:2];
        dq3 = (upper_Bound-2)/N;
        q_back3 = [2:dq3:upper_Bound];
        dq1_simp = ones(size(q_back1,2),1)';%%% the formula of trapiech
        dq1_simp(1) = 1/2;
        dq1_simp(size(q_back1,2)) = 1/2;
        dq2_simp = ones(size(q_back2,2),1)';%%% the formula of trapiech
        dq2_simp(1) = 1/2;
        dq2_simp(size(q_back2,2)) = 1/2;
        dq3_simp = ones(size(q_back3,2),1)';%%% the formula of trapiech
        dq3_simp(1) = 1/2;
        dq3_simp(size(q_back3,2)) = 1/2;

        q_back = [q_back1,q_back2, q_back3];
        dq_simp = [dq * dq1_simp, dq2 * dq2_simp, dq3 * dq3_simp];
    else
        dq = (2-q_start)/N;
        q_back1 = [[q_start:dq:1-dq],[1+dq:dq:2]];
        dq2 = (upper_Bound - 2)/N;
        q_back2 = [2:dq2:upper_Bound];
        dq1_simp = ones(size(q_back1,2),1)';%%% the formula of trapiech
        dq1_simp(1) = 1/2;
        dq1_simp(size(q_back1,2)) = 1/2;
        dq2_simp = ones(size(q_back2,2),1)';%%% the formula of trapiech
        dq2_simp(1) = 1/2;
        dq2_simp(size(q_back2,2)) = 1/2;

        q_back = [q_back1,q_back2];
        dq_simp = [dq * dq1_simp, dq2 * dq2_simp];
    end

    
    case 'simpson'
    
%%%%%%%%% Simpson formula start
if(q_0/1i ~= imag(q_0) && q_0 > q_start)
    dq1 = (q_0-q_start)/N;
    q_back1 = [q_start:dq1:q_0-dq1];
    dq2 = (2-q_0-dq1)/N;
    q_back2 = [q_0+dq1:dq2:2];
    if(rem(size(q_back2,2),2))
        q_back2 = [q_back2, 2+dq2];
    end
    dq3 = (upper_Bound-q_back2(end))/N;
    q_back3 = [q_back2(end):dq3:upper_Bound];
    if(rem(size(q_back3,2),2))
        q_back3 = [q_back3, upper_Bound+dq3];
    end
    
    dq1_simp = ones(size(q_back1,2),1)';
    dq1_simp(1:2:end) = 2;
    dq1_simp(2:2:end) = 4;
    dq1_simp(1) = 1;
    dq1_simp(size(q_back1,2)) = 1;
    dq2_simp = ones(size(q_back2,2),1)';
    dq2_simp(1:2:end) = 2;
    dq2_simp(2:2:end) = 4;
    dq2_simp(1) = 1;
    dq2_simp(size(q_back2,2)) = 1;
    dq3_simp = ones(size(q_back3,2),1)';
    dq3_simp(1:2:end) = 2;
    dq3_simp(2:2:end) = 4;
    dq3_simp(1) = 1;
    dq3_simp(size(q_back3,2)) = 1;

    dq1_simp = dq1_simp.* dq1;
    dq2_simp = dq2_simp.* dq2;
    dq3_simp = dq3_simp.* dq3;
    dq = 1/3;

    q_back = [q_back1,q_back2, q_back3];
    dq_simp = dq * [dq1_simp, dq2_simp, dq3_simp];
 else
%     dq1 = (2-q_start)/N;
% %     q_back1 = [q_start:dq1:2];
%     q_back1 = [[q_start:dq1:1-dq1],[1+dq1:dq1:2]];
%     dq2 = (upper_Bound-2)/N;
%     q_back2 = [2:dq2:upper_Bound];
%     if(rem(size(q_back2,2),2))
%         q_back2 = [q_back2, upper_Bound+dq2];
%     end
%     
%     dq1_simp = ones(size(q_back1,2),1)';
%     dq1_simp(1:2:end) = 2;
%     dq1_simp(2:2:end) = 4;
%     dq1_simp(1) = 1;
%     dq1_simp(size(q_back1,2)) = 1;
%     dq2_simp = ones(size(q_back2,2),1)';
%     dq2_simp(1:2:end) = 2;
%     dq2_simp(2:2:end) = 4;
%     dq2_simp(1) = 1;
%     dq2_simp(size(q_back2,2)) = 1;
% 
%     dq1_simp = dq1_simp.* dq1;
%     dq2_simp = dq2_simp.* dq2;
%     dq = 1/3;
% 
%     q_back = [q_back1,q_back2];
%     dq_simp = dq * [dq1_simp, dq2_simp];
   
    
    dq1 = (1-1e-8 -q_start)/N;
    q_back1 = [q_start:dq1:1-1e-8];
    dq2 = (2-(1+1e-8)-dq1)/N;
    q_back2 = [1+1e-8:dq2:2];
    if(rem(size(q_back2,2),2))
        q_back2 = [q_back2, 2+dq2];
    end
    dq3 = (upper_Bound-q_back2(end) - dq2)/N;
    q_back3 = [q_back2(end)+ dq2:dq3:upper_Bound];
    if(rem(size(q_back3,2),2))
        q_back3 = [q_back3, upper_Bound+dq3];
    end
    
    dq1_simp = ones(size(q_back1,2),1)';
    dq1_simp(1:2:end) = 2;
    dq1_simp(2:2:end) = 4;
    dq1_simp(1) = 1;
    dq1_simp(size(q_back1,2)) = 1;
    dq2_simp = ones(size(q_back2,2),1)';
    dq2_simp(1:2:end) = 2;
    dq2_simp(2:2:end) = 4;
    dq2_simp(1) = 1;
    dq2_simp(size(q_back2,2)) = 1;
    dq3_simp = ones(size(q_back3,2),1)';
    dq3_simp(1:2:end) = 2;
    dq3_simp(2:2:end) = 4;
    dq3_simp(1) = 1;
    dq3_simp(size(q_back3,2)) = 1;

    dq1_simp = dq1_simp.* dq1;
    dq2_simp = dq2_simp.* dq2;
    dq3_simp = dq3_simp.* dq3;
    dq = 1/3;

    q_back = [q_back1,q_back2, q_back3];
    dq_simp = dq * [dq1_simp, dq2_simp, dq3_simp];    

        %%%%%%%%%% Simpson formula end
end
    case'difficult_simpson'
     dq = 0.01;
    N1 = 50;
    N2 = 25;
    N3 = 400;
    
    dq1 = (q_0-dq-1e-7)/N1;
    dq2 = (1-1e-7-q_0-dq)/N2;
    dq3 = (upper_Bound-1-1e-7)/N3;
    
    q_back1 = [1e-7:dq1:q_0-dq];
    q_back2 = [q_0+dq:dq2:1-1e-7];
    q_back3 = [1+1e-7:dq3:upper_Bound];
    dq1_simp = ones(size(q_back1,2),1)';
    dq1_simp(1:2:end) = 2;
    dq1_simp(2:2:end) = 4;
    dq1_simp(1) = 1;
    dq1_simp(size(q_back1,2)) = 1;
    dq2_simp = ones(size(q_back2,2),1)';
    dq2_simp(1:2:end) = 2;
    dq2_simp(2:2:end) = 4;
    dq2_simp(1) = 1;
    dq2_simp(size(q_back2,2)) = 1;
    dq3_simp = ones(size(q_back3,2),1)';
    dq3_simp(1:2:end) = 2;
    dq3_simp(2:2:end) = 4;
    dq3_simp(1) = 1;
    dq3_simp(size(q_back3,2)) = 1;

    dq1_simp = dq1_simp.* dq1;
    dq2_simp = dq2_simp.* dq2;
    dq3_simp = dq3_simp.* dq3;
    dq = 1/3;

    q_back = [q_back1,q_back2,q_back3];
    dq_simp = dq * [dq1_simp, dq2_simp, dq3_simp];
    %%%%%%%%%% Simpson formula end
end
            
    p = sqrt(1 - q_back.^2);
    p = real(p) - 1i * abs(imag(p));
