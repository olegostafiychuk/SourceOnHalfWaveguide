function [dq_simp, q_back, p] = quadratureMethod_forIntegralEqs_multipleIntervals(method, N, N_upper, q_n, q_start, upper_Bound)


switch(method)
    case'difficult_simpsonMarkovGildendurgExpWithCollisions'
    
    dq = 1e-4;
    
    if(isempty(q_n))
        N1 = 200;
        N2 = 800;
        N3 = 800;
        N4 = 1600;
    
        middleBound = 30000;
        middleBound2 = 300000;
    
        dq1 = (1-dq -q_start)/N1;
        dq2 = (middleBound-(1+dq))/N2;
        dq3 = (middleBound2-(middleBound+dq))/N3;
        dq4 = (upper_Bound-(middleBound2+dq2))/N4;
    
        q_back1 = [q_start:dq1:1-dq];
        q_back2 = [1+dq:dq2:middleBound];
        q_back3 = [middleBound+dq:dq3:middleBound2];
        q_back4 = [middleBound2+dq2:dq4:upper_Bound];
        %     if(rem(size(q_back2,2),2))
        %         q_back2 = [q_back2, 2+dq2];
        %     end
        %     if(rem(size(q_back4,2),2))
        %         q_back4 = [q_back4, upper_Bound+dq4];
        %     end

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
        dq4_simp = ones(size(q_back4,2),1)';
        dq4_simp(1:2:end) = 2;
        dq4_simp(2:2:end) = 4;
        dq4_simp(1) = 1;
        dq4_simp(size(q_back4,2)) = 1;
        
        dq1_simp = dq1_simp.* dq1;
        dq2_simp = dq2_simp.* dq2;
        dq3_simp = dq3_simp.* dq3;
        dq4_simp = dq4_simp.* dq4;
        dq = 1/3;
        
        q_back = [q_back1, q_back2, q_back3, q_back4];
        dq_simp = dq * [dq1_simp, dq2_simp, dq3_simp, dq4_simp];
        %%%%%%%%%% Simpson formula end
    else
        N1 = 400;
%         dq1 = (1-dq -q_start) / N1;
%         dq2 = (q_n(1)-(1+dq)) / N1;
%         dq1 = (1-dq -q_start) / N_upper;
        dq1 = (1-dq -q_start) / N1;
        dq2 = (q_n(1)-(1+dq)) / N;
        q_back1 = [q_start:dq1:1-dq];
        q_back2 = [1+dq:dq2:q_n(1)];
        
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
        
        q_back = [q_back1, q_back2];
        dqq = 1/3;
        dq1_simp = dq1_simp.* dq1;
        dq2_simp = dq2_simp.* dq2;
        dq_simp = dqq * [dq1_simp, dq2_simp];
        
        
        dq2 = dq;
        for iq = 2:size(q_n)-1
            dq3 = (q_n(iq) - dq2 - q_n(iq-1) - dq2) / N;
            q_back3 = [q_n(iq-1) + dq2:dq3:q_n(iq) - dq2];
            
            dq3_simp = ones(size(q_back3,2),1)';
            dq3_simp(1:2:end) = 2;
            dq3_simp(2:2:end) = 4;
            dq3_simp(1) = 1;
            dq3_simp(size(q_back3,2)) = 1;
            
            dq3_simp = dq3_simp.* dq3;
            q_back = [q_back, q_back3];
            dq_simp =  [dq_simp, dqq * dq3_simp];
%             dq2 = dq3; 
        end
        
        dq4 = (q_n(end) - dq2 - (q_n(end-1) + dq2)) / N_upper;
        q_back4 = [q_n(end-1) + dq2:dq4:q_n(end) - dq2];
        dq4_simp = ones(size(q_back4,2),1)';
        dq4_simp(1:2:end) = 2;
        dq4_simp(2:2:end) = 4;
        dq4_simp(1) = 1;
        dq4_simp(size(q_back4,2)) = 1;
        dq4_simp = dq4_simp.* dq4;
        
        q_back = [q_back, q_back4];
        dq_simp =  [dq_simp, dqq * dq4_simp];
        
        dq3 = dq;        
        dq4 = (upper_Bound - (q_n(end) + dq3)) / N_upper;
        q_back4 = [q_n(end) + dq3:dq4:upper_Bound];
        dq4_simp = ones(size(q_back4,2),1)';
        dq4_simp(1:2:end) = 2;
        dq4_simp(2:2:end) = 4;
        dq4_simp(1) = 1;
        dq4_simp(size(q_back4,2)) = 1;
        dq4_simp = dq4_simp.* dq4;
        
        q_back = [q_back, q_back4];
        dq_simp =  [dq_simp, dqq * dq4_simp];
    end
        %%%%%%%%%% Simpson formula end
    case'simple_simpson_ExpWithCollisions'
    
        dq1 = (upper_Bound -q_start)/ N_upper;    
        q_back1 = [q_start:dq1:upper_Bound];

        dq1_simp = ones(size(q_back1,2),1)';
        dq1_simp(1:2:end) = 2;
        dq1_simp(2:2:end) = 4;
        dq1_simp(1) = 1;
        dq1_simp(size(q_back1,2)) = 1;        
        dq1_simp = dq1_simp.* dq1;
        dq = 1/3;

        q_back = [q_back1];
        dq_simp = dq * [dq1_simp];
end
            
    p = sqrt(1 - q_back.^2);
    p = p.* (2*(imag(p) <= 0)-1);
