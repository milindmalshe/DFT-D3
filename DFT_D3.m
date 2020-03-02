%function [E_disp, S] = Milind_DFT_D3(S)
function [E_disp, S] = DFT_D3()

% Number of elements for which DFT-D3 was parameterized
NumElements = 94;

% Coordination number - number of reference systems considered for 2-D
% interpolation
CoordinationNum_ref_max = 5;

Damping = 0;

numAtoms = 2;

% Steepness parameters alpha
alpha6 = 14;
alpha8 = alpha6 + 2; %16;
alpha = [alpha6, alpha8];

k1 = 16.0;
k2 = 4.0/3.0;
k3 = -4.0;

 
[CN_A_Num, ElementType] = CoordinationStaeReference_CN();
[CN_A] = CoordinationNumbers_CN();
[r2_r4] = r_cut();

% R_Acov = Sclaed covalent radius of atom A
[R_cov] = rcov();

if (Damping ==0)
% Stefan Grimme Zero Damping
    Damp=0;
elseif (Damping == 1)
% Use Becke-Johnson Damping
    Damp = 1; 
end

[RS6,S18,RS18,S6] = D3_DampingParam(Damp);




% CN_A = Fractional Coordination Number 
CN_A = [];

[c6_AB] = c6();
num_atoms=1;

r = 1;

E_disp = 0;
E2 = 0;
E3 = 0;

% Global (density functional dependent) scaling factors Sn
% Sn are adjusted only for n>6 to ensure aymptotic exactness: S6 = 1;
% Sn = [1,0,8];
% S6 = Sn(1,1);
% S7 = Sn(1,2);
% S8 = Sn(1,3);
s6 = 1;
s7 = 0;
s8 = 0.5;
scaling_fact = [s6, s7, s8];

% n-th order dispersion coefficients
% higher-order dispersion coefficients are computed recursively
C6 = 1;
C8 = 1;
disp_coeff = [C6, C8];


%Order-dependent scaling factor of the cut-off radii R_cutOff
Sr8 = 1;
Sr7 = 0;
Sr6 = 1; % $$$$ CHECK VALUE OF S_r6


% Order dependenet scaling factor of cutoff radii Srn

Srn = [Sr6, Sr7, Sr8];







% Compute E2
% E2 = -0.5;
r = [1, 1; 1, 1];

for i= 1:numAtoms
    for j= 1:numAtoms
        if(j ~= i)
            % Compute damping function fd,n
            fd6= 1+ 6.*((r(i,j)./(Sr6.*r(i,j))).^(-alpha6));
            fd6= 1./fd6;

            fd8= 1+ 6.*((r(i,j)./(Sr8.*r(i,j))).^(-alpha8));
            fd8= 1./fd8;
            
            

            E26 = (s6.*(C6./r(i,j).^6).*fd6);
            E28 = (s8.*(C8./r(i,j).^8).*fd8); 
        
            E2 = E26+E28;
            
            
        end
        
    end
    
end


[CN_A] = CN(numAtoms, k1, k2, R_cov, r);


E_disp = E2 + E3;








n_disp_coeff= [6,8];
num_disp_coeff = size(n_disp_coeff,2);

% % damping_func = 1/(1+6*(r/(SSR_cutOff) )^(-alpha) );

%E2
