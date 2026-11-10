function [A, C, C6] = murnaghan_A(lam, mu, l, m, n, sig, eps)
% MURNAGHAN_A  参考构型上的一阶声弹有效模量张量（自然速度用）
%   A_ijpn = C_ijpn + σ°_nj δ_ip + C_ijpnqr ε°_qr + ε°_ik C_kjpn + ε°_pm C_ijmn
% 输入：λ, μ 二阶；l, m, n Murnaghan 三阶（Hughes–Kelly 定义，拉为正）；
%       sig 3×3 预应力（≈ 二阶 PK 应力，一阶下等于 Cauchy 应力），eps 3×3 预应变。
% 输出：A(i,j,p,n)：增量第一 PK 应力 δP_ij = A_ijpn ∂u_p/∂X_n；
%       C 二阶张量；C6 三阶张量（W = ½C:E:E + ⅙C6:E:E:E）。
% 与 Landau 常数的关系：A_L = n, B_L = m − n/2, C_L = l − m + n/2。
% A 有大对称 (ij)<->(pn)，没有小对称（i<->j 不对称），所以要用位移梯度而不是应变。

d = eye(3);
AL = n;  BL = m - n/2;  CL = l - m + n/2;
I4 = @(i,j,k,l_) d(i,k)*d(j,l_) + d(i,l_)*d(j,k);

C = zeros(3,3,3,3);
for i = 1:3, for j = 1:3, for k = 1:3, for l_ = 1:3
    C(i,j,k,l_) = lam*d(i,j)*d(k,l_) + mu*I4(i,j,k,l_);
end, end, end, end

C6 = zeros(3,3,3,3,3,3);
for i = 1:3, for j = 1:3, for k = 1:3, for l_ = 1:3, for m_ = 1:3, for n_ = 1:3
    C6(i,j,k,l_,m_,n_) = 2*CL*d(i,j)*d(k,l_)*d(m_,n_) ...
        + BL*( d(i,j)*I4(k,l_,m_,n_) + d(k,l_)*I4(i,j,m_,n_) + d(m_,n_)*I4(i,j,k,l_) ) ...
        + (AL/4)*( d(i,k)*I4(j,l_,m_,n_) + d(i,l_)*I4(j,k,m_,n_) ...
                 + d(j,k)*I4(i,l_,m_,n_) + d(j,l_)*I4(i,k,m_,n_) );
end, end, end, end, end, end

A = zeros(3,3,3,3);
for i = 1:3, for j = 1:3, for p = 1:3, for q = 1:3      % q 即公式里的 n
    s = C(i,j,p,q) + sig(q,j)*d(i,p);
    for a = 1:3, for b = 1:3
        s = s + C6(i,j,p,q,a,b)*eps(a,b);
    end, end
    for k = 1:3
        s = s + eps(i,k)*C(k,j,p,q) + eps(p,k)*C(i,j,k,q);
    end
    A(i,j,p,q) = s;
end, end, end, end
end
