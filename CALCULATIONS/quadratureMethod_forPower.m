function [dr_simp, r1] = quadratureMethod_forPower(method, N, r_start, upper_Bound)


switch(method)
    case 'simpson'
    
%%%%%%%%% Simpson formula start
    
    dr1 = (upper_Bound - r_start)/N;
    r1 = [r_start:dr1:upper_Bound];
   
    dr1_simp = ones(size(r1,2),1)';
    dr1_simp(1:2:end) = 2;
    dr1_simp(2:2:end) = 4;
    dr1_simp(1) = 1;
    dr1_simp(size(r1,2)) = 1;
    dr1_simp = dr1_simp.* dr1;
    dr = 1/3;
    dr_simp = dr * dr1_simp;
 
end
