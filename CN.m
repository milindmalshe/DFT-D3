function [CN_A] = CN(numAtoms, k1, k2, R_cov, r)

% Calculate the coordination number of each atom in the system

for i=1:numAtoms
    for j=1:numAtoms
        if (j ~=i )
           R_cov_AB = R_cov(i)+R_cov(j);
           
           CN_exp= exp(-k1.*(k2.*R_cov_AB./r(i,j) - 1));
           CN_A(i)= 1./(1+CN_exp);
           
        end
        
    end
end
