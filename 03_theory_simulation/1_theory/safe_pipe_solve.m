function modes = safe_pipe_solve(S, omega, D, rho_s, Ks, rho_f, cmin, cmax)
% SAFE_PIPE_SOLVE  给定 ω，解轴对称充液管的实波数模态并做后处理
%   D     : 5×5（均匀）或 5×5×nGs（每个壁 Gauss 点一个，预应力直接法用）
%   rho_s : 壁密度；Ks, rho_f : 液体绝热体积模量与密度
%   cmin, cmax : 只保留相速度在此范围内、k 为实正的模态
% 二次特征问题 (A0 + i k A1 + k^2 A2) q = 0：
%   壁：K0 + i k K1 + k^2 K2 − ω^2 M，h = (B0 + i k B1) q_e，
%       B0 行：H_rr = dN u_r，H_θθ = N u_r / r，H_zr = dN u_z；B1 行：H_zz = N u_z，H_rz = N u_r
%   液：(1/ρ_f)(p'q' + k^2 p q) − (ω^2/K_s) p q；界面：壁受 +p e_r，液受 −ω^2 r_i u_r
% 输出 modes(m)：k, cp, cg, Vs, Ts, Vf, Tf, Pw, eta_f(液体能量占比), q(整体向量),
%   以及 Gauss 点场：h (5×nGs), u (2×nGs), p (nGf×1), dp (nGf×1)
n = S.ndof;  Ns = S.Ns;
K0 = zeros(n); K1 = zeros(n); K2 = zeros(n); M = zeros(n);
uniformD = (ndims(D) == 2);
%% 壁
g = S.s;
for i = 1:numel(g.r)
    r = g.r(i); w = g.w(i)*r; N = g.N(i,:); dN = g.dN(i,:);
    if uniformD, Di = D; else, Di = D(:,:,i); end
    B0 = zeros(5,6); B1 = zeros(5,6);
    B0(1,1:2:5) = dN;  B0(2,1:2:5) = N/r;  B0(5,2:2:6) = dN;
    B1(3,2:2:6) = N;   B1(4,1:2:5) = N;
    Nu = zeros(2,6); Nu(1,1:2:5) = N; Nu(2,2:2:6) = N;
    dof = reshape([2*g.nodes(i,:)-1; 2*g.nodes(i,:)], 1, 6);
    K0(dof,dof) = K0(dof,dof) + w*(B0'*Di*B0);
    K1(dof,dof) = K1(dof,dof) + w*(B0'*Di*B1 - B1'*Di*B0);
    K2(dof,dof) = K2(dof,dof) + w*(B1'*Di*B1);
    M(dof,dof)  = M(dof,dof)  + w*rho_s*(Nu'*Nu);
end
%% 液
g = S.f;
for i = 1:numel(g.r)
    r = g.r(i); w = g.w(i)*r; N = g.N(i,:); dN = g.dN(i,:);
    dof = 2*Ns + g.nodes(i,:);
    K0(dof,dof) = K0(dof,dof) + w*(dN'*dN)/rho_f;
    K2(dof,dof) = K2(dof,dof) + w*(N'*N)/rho_f;
    M(dof,dof)  = M(dof,dof)  + w*(N'*N)/Ks;
end
%% 耦合与缩放
% 未缩放的耦合是非对称的（壁行 −r_i p，液行 −ω^2 r_i u_r），且壁块与液块量级差 ~1e13，polyeig 失准。
% 令 p = β ω^2 p~，液体方程乘 β：液块 → β^2 ω^2 K_f，两个耦合项都成 −β ω^2 r_i，系统变 Hermitian。
% β 取使液块与壁块 Frobenius 范数相当。
ri = S.geom.r_i;
is = 1:2*Ns;  ifl = 2*Ns + (1:S.Nf);
A0 = K0 - omega^2*M;
beta = sqrt(norm(A0(is,is), 'fro')/norm(A0(ifl,ifl), 'fro'))/omega;
sc = ones(n,1);  sc(ifl) = beta*omega;            % 对称缩放向量：A_scaled = diag(sc) A diag(sc)
A0 = (sc*sc').*A0;  K1s = (sc*sc').*K1;  K2s = (sc*sc').*K2;
A0(S.i_ur_i, S.i_p_i) = A0(S.i_ur_i, S.i_p_i) - beta*omega^2*ri;
A0(S.i_p_i, S.i_ur_i) = A0(S.i_p_i, S.i_ur_i) - beta*omega^2*ri;
%% 解（再把 k 无量纲化：k = κ λ，κ = ω/c_T）
kap = omega/sqrt(D(4,4,1)/rho_s);
s0  = norm(A0, 'fro');
[X, e] = polyeig(A0/s0, 1i*K1s*kap/s0, K2s*kap^2/s0);
e = e*kap;
sel = abs(imag(e)) < 1e-6*abs(real(e)) & real(e) > 0;
k = real(e(sel));  X = X(:, sel);
X(ifl, :) = X(ifl, :)*beta*omega^2;               % 还原 p = β ω^2 p~
cp = omega./k;
keep = cp >= cmin & cp <= cmax;
k = k(keep); X = X(:, keep); cp = cp(keep);
[cp, ord] = sort(cp, 'descend'); k = k(ord); X = X(:, ord);
%% 后处理
modes = struct([]);
for m = 1:numel(k)
    q = X(:, m);  km = k(m);
    md.k = km; md.cp = cp(m); md.q = q;
    % 壁场
    gs = S.s; nG = numel(gs.r);
    h = zeros(5, nG); u = zeros(2, nG);
    Vs = 0; Ts = 0; Ps = 0;
    for i = 1:nG
        r = gs.r(i); w = gs.w(i)*r; N = gs.N(i,:); dN = gs.dN(i,:);
        if uniformD, Di = D; else, Di = D(:,:,i); end
        dof = reshape([2*gs.nodes(i,:)-1; 2*gs.nodes(i,:)], 1, 6);
        qe = q(dof);
        B0 = zeros(5,6); B1 = zeros(5,6);
        B0(1,1:2:5) = dN;  B0(2,1:2:5) = N/r;  B0(5,2:2:6) = dN;
        B1(3,2:2:6) = N;   B1(4,1:2:5) = N;
        hi = (B0 + 1i*km*B1)*qe;
        ui = [N*qe(1:2:5); N*qe(2:2:6)];
        h(:,i) = hi; u(:,i) = ui;
        Vs = Vs + 0.25*w*real(hi'*Di*hi);
        Ts = Ts + 0.25*w*rho_s*omega^2*(ui'*ui);
        s = Di*hi;                                  % 增量应力 (rr,θθ,zz,rz,zr)
        Ps = Ps + 0.5*omega*w*imag(conj(ui(1))*s(4) + conj(ui(2))*s(3));
    end
    % 液场
    gf = S.f; nGf = numel(gf.r);
    p = zeros(nGf,1); dp = zeros(nGf,1);
    Vf = 0; Tf = 0; Pf = 0;
    for i = 1:nGf
        r = gf.r(i); w = gf.w(i)*r; N = gf.N(i,:); dN = gf.dN(i,:);
        qe = q(2*Ns + gf.nodes(i,:));
        pi_ = N*qe; dpi = dN*qe;
        p(i) = pi_; dp(i) = dpi;
        Vf = Vf + 0.25*w*abs(pi_)^2/Ks;
        Tf = Tf + 0.25*w*(abs(dpi)^2 + km^2*abs(pi_)^2)/(omega^2*rho_f);
        Pf = Pf + 0.5*w*km*abs(pi_)^2/(omega*rho_f);
    end
    md.h = h; md.u = u; md.p = p; md.dp = dp;
    md.Vs = Vs; md.Ts = Ts; md.Vf = Vf; md.Tf = Tf;
    md.E = Vs + Ts + Vf + Tf;
    md.Pw = Ps + Pf;
    md.cg = md.Pw/md.E;
    md.eta_f = (Vf + Tf)/md.E;
    md.equip = (Ts + Tf)/(Vs + Vf) - 1;            % 均分检验，应 ≈ 0
    modes = [modes, md]; %#ok<AGROW>
end
end
