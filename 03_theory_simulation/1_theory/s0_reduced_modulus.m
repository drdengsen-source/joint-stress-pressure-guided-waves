function rhoV2 = s0_reduced_modulus(A)
% S0_REDUCED_MODULUS  薄板极限下 S0 板波的 ρ0·V_nat^2
%   传播沿 X1，横向 X2 平面应变（∂/∂X2 = 0），面外 X3 两表面增量牵引为零。
%   ρ0 V^2 = A_1111 − Σ_{p,k} A_11p3 (Q^-1)_pk A_k311,   Q_ip = A_i3p3
%   无预应力各向同性时退回 E/(1−ν^2)。
Q = zeros(3,3);  R = zeros(3,1);  S = zeros(1,3);
for i = 1:3
    for p = 1:3, Q(i,p) = A(i,3,p,3); end
    R(i) = A(i,3,1,1);
    S(i) = A(1,1,i,3);
end
rhoV2 = A(1,1,1,1) - S*(Q\R);
end
