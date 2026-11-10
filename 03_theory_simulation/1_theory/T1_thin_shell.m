%% T1_thin_shell.m  薄壳闭式极限
% 内容：
%   1. 壁模态：S0 板波在双轴面内应力 (σz, σθ) 下相速度变化的闭式 → K_z, K_θ（自然／真实）
%      附 σr 交叉项估计，与"每 MPa 内压"的交叉灵敏度 K_p
%   2. 体波退化：用同一张量算 Hughes–Kelly 五个体波公式，证明 H–K 公式是真实速度
%   3. 液模态：Korteweg + c_f(p,T) → (1/c_K)dc_K/dp、dc_K/dT，κ(0)
%   4. 约定三：几何更新对两族灵敏度的贡献
%   5. 汇总表（延时灵敏度 ns/MPa、ns/K）
% 约定：拉为正；Murnaghan (l,m,n) 按 Hughes–Kelly；自然速度 c_nat 由参考构型算，
%       真实速度 c_true = c_nat (1+εz)；EMAT 固定标距 W_E 的延时 t = W_E/c_true。
% 用法：直接运行。换牌号／液体改 params_pipe.m。结果同时写到 T1_结果.txt。
clear; clc;
here = fileparts(mfilename('fullpath'));
P = params_pipe('7075', 'oil');
W = P.wall;  F = P.fluid;  G = P.geom;
fid = fopen(fullfile(here, 'T1_结果.txt'), 'w', 'n', 'UTF-8');
pr(fid, 'T1 薄壳闭式极限  %s\n', datestr(now, 'yyyy-mm-dd HH:MM'));
pr(fid, '管壁：%s\n', W.name);
pr(fid, '  λ=%.3f  μ=%.3f  l=%.1f  m=%.1f  n=%.1f GPa   ρ=%g kg/m^3\n', ...
    W.lam/1e9, W.mu/1e9, W.l/1e9, W.m/1e9, W.n/1e9, W.rho);
pr(fid, '  E=%.2f GPa  ν=%.4f  cL=%.0f  cT=%.0f  c_plate=%.0f  c_bar=%.0f m/s\n', ...
    W.E/1e9, W.nu, W.cL, W.cT, W.c_plate, W.c_bar);
pr(fid, '液体：%s  c_f=%g m/s  ρ_f=%g  B/A=%g  K_s=%.3f GPa\n', F.name, F.c, F.rho, F.BA, F.Ks/1e9);
pr(fid, '几何：r_i=%.1f r_o=%.1f t=%.1f mm  W_E=%g mm  A_wall=%.2f  A_bore=%.2f mm^2\n', ...
    G.r_i*1e3, G.r_o*1e3, G.t*1e3, G.W_E*1e3, G.A_wall*1e6, G.A_bore*1e6);
pr(fid, '内压换算：σθ/p=%.2f  σz(闭端)/p=%.3f  σr(均值)/p=%.2f\n\n', ...
    G.sig_theta_per_p, G.sig_z_closed_per_p, G.sig_r_per_p);

MPa = 1e6;
hooke = @(sig) ((1 + W.nu)*sig - W.nu*trace(sig)*eye(3))/W.E;   % 各向同性预应变

%% 1. 壁模态 S0：三个方向单位应力的自然系数
% 坐标：1 = z（传播），2 = θ（横向面内），3 = r（面外）
A0 = murnaghan_A(W.lam, W.mu, W.l, W.m, W.n, zeros(3), zeros(3));
rhoV2_0 = s0_reduced_modulus(A0);
c_S0 = sqrt(rhoV2_0/W.rho);
pr(fid, '== 1. 壁模态（薄板 S0 极限） ==\n');
pr(fid, 'ρ0 V0^2 = %.3f GPa  (E/(1-ν^2) = %.3f GPa)   V0 = %.1f m/s\n', ...
    rhoV2_0/1e9, W.E/(1 - W.nu^2)/1e9, c_S0);

K_nat = zeros(1, 3);
for dir = 1:3
    sig = zeros(3);  sig(dir, dir) = 1*MPa;
    A = murnaghan_A(W.lam, W.mu, W.l, W.m, W.n, sig, hooke(sig));
    K_nat(dir) = (s0_reduced_modulus(A) - rhoV2_0)/(2*rhoV2_0);    % 每 MPa
end
% 线性度检查（100 MPa）
sig = zeros(3); sig(1,1) = 100*MPa;
A = murnaghan_A(W.lam, W.mu, W.l, W.m, W.n, sig, hooke(sig));
K100 = (sqrt(s0_reduced_modulus(A)/rhoV2_0) - 1)/100;
pr(fid, '自然系数 (1/V)dV/dσ  [1/MPa]:  K_z=%.4e  K_θ=%.4e  K_r(估计)=%.4e\n', K_nat);
pr(fid, '  线性度：100 MPa 单轴按 sqrt 精确算 %.4e，与一阶 %.4e 差 %.2f%%\n', ...
    K100, K_nat(1), (K100/K_nat(1) - 1)*100);

% 真实速度：c_true = c_nat (1+εz)，εz = [σz − ν(σθ+σr)]/E
K_true = K_nat + [1, -W.nu, -W.nu]/(W.E/MPa);
pr(fid, '真实系数 [1/MPa]:            K_z=%.4e  K_θ=%.4e  K_r(估计)=%.4e\n', K_true);
pr(fid, '  εz 修正 1/E = %.3e /MPa，占 K_z 的 %.1f%%\n', 1/(W.E/MPa), 100/(W.E/MPa)/abs(K_nat(1)));

% 每 MPa 内压的交叉项（σz 总量单独作未知量；这里只含 σθ、σr）
Kp_nat  = G.sig_theta_per_p*K_nat(2)  + G.sig_r_per_p*K_nat(3);
Kp_true = G.sig_theta_per_p*K_true(2) + G.sig_r_per_p*K_true(3);
Kp_true_noR = G.sig_theta_per_p*K_true(2);
pr(fid, '内压交叉项 K_p = (r_i/t)K_θ + (−1/2)K_r [1/MPa]:  自然 %.4e   真实 %.4e   (不计 σr：%.4e)\n', ...
    Kp_nat, Kp_true, Kp_true_noR);
pr(fid, '  单模态法忽略交叉项的假应力：%.3f MPa / MPa 内压（真实速度）\n', Kp_true/K_true(1));

% 延时灵敏度
t0_w = G.W_E/c_S0;
dt_dsz = -t0_w*K_true(1);   dt_dp_w = -t0_w*Kp_true;
pr(fid, '延时 t0 = %.2f µs：∂t/∂σz = %+.3f ns/MPa，∂t/∂p(交叉) = %+.3f ns/MPa\n\n', ...
    t0_w*1e6, dt_dsz*1e9, dt_dp_w*1e9);

%% 2. 体波退化：Hughes–Kelly 校核（单轴拉伸 T，传播 N，极化 a）
pr(fid, '== 2. 体波退化：同一张量 vs Hughes–Kelly 公式（T = 100 MPa 单轴，沿 1）==\n');
T = 100*MPa;  sig = zeros(3); sig(1,1) = T;  eps = hooke(sig);
A = murnaghan_A(W.lam, W.mu, W.l, W.m, W.n, sig, eps);
lam = W.lam; mu = W.mu; l = W.l; m = W.m; n = W.n; K3 = 3*W.K0;
HK = struct();   % 真实速度形式（原文压为正的 P 换成 −T）
HK.v11 = lam + 2*mu + T/K3*(2*l + lam + (lam + mu)/mu*(4*m + 4*lam + 10*mu));
HK.v22 = lam + 2*mu + T/K3*(2*l - 2*lam/mu*(m + lam + 2*mu));
HK.v12 = mu + T/K3*(m + lam*n/(4*mu) + 4*lam + 4*mu);
HK.v21 = mu + T/K3*(m + lam*n/(4*mu) + lam + 2*mu);
HK.v23 = mu + T/K3*(m - (lam + mu)/(2*mu)*n - 2*lam);
cases = {'v11', 1, 1; 'v22', 2, 2; 'v12', 1, 2; 'v21', 2, 1; 'v23', 2, 3};   % {名, 传播 N, 极化 a}
pr(fid, '%-6s %14s %14s %14s %10s\n', '波', 'ρ0V^2 自然', 'ρ0v^2 真实', 'H–K 公式', '真实/HK−1');
for ic = 1:size(cases, 1)
    N = zeros(3,1); N(cases{ic,2}) = 1;  a = cases{ic,3};
    Gam = zeros(3);
    for i = 1:3, for p = 1:3, for j = 1:3, for q = 1:3
        Gam(i,p) = Gam(i,p) + A(i,j,p,q)*N(j)*N(q);
    end, end, end, end
    rhoV2_nat = Gam(a,a);                                       % 主轴对齐，Γ 对角
    epsN = eps(cases{ic,2}, cases{ic,2});
    rhoV2_true = rhoV2_nat*(1 + epsN)^2;
    hk = HK.(cases{ic,1});
    pr(fid, '%-6s %14.4f %14.4f %14.4f %10.2e\n', cases{ic,1}, rhoV2_nat/1e9, rhoV2_true/1e9, hk/1e9, rhoV2_true/hk - 1);
end
pr(fid, '结论：真实速度形式与 H–K 五式逐一吻合（残差 ≤ 4e−5 是 (1+ε)^2 与 1+2ε 之差乘一阶变化，二阶量），H–K 公式是真实速度。\n');
% 体波相对灵敏度（每 MPa）
k11_true = (sqrt(HK.v11/(lam + 2*mu)) - 1)/100;   k22_true = (sqrt(HK.v22/(lam + 2*mu)) - 1)/100;
pr(fid, '体波纵波真实系数：沿应力 %.3e /MPa，垂直应力 %.3e /MPa\n\n', k11_true, k22_true);

%% 3. 液模态：Korteweg + c_f(p, T)
pr(fid, '== 3. 液模态（Korteweg 长波极限） ==\n');
Om_thin = 2*G.r_i*F.Ks/(W.E*G.t);                                   % 薄环：εθ(r_i)/p = r_i/(E t)
eps_hat = ((G.r_o^2 + G.r_i^2)/(G.r_o^2 - G.r_i^2) + W.nu)/W.E;     % 厚壁 Lamé（轴向自由）：εθ(r_i)/p
Om = 2*F.Ks*eps_hat;                                                % Ω = 2 K_s εθ(r_i)/p
c_K = F.c/sqrt(1 + Om);
pr(fid, 'Ω 薄环 2 r_i K_s/(E t) = %.4f → c_K = %.1f；Ω 厚壁 2 K_s[(r_o²+r_i²)/(r_o²−r_i²)+ν]/E = %.4f → c_K = %.1f m/s（c_f = %g）\n', ...
    Om_thin, F.c/sqrt(1 + Om_thin), Om, c_K, F.c);
pr(fid, '  （SAFE 精确解 5 kHz 给 1262 m/s，取厚壁公式；以下都用厚壁 Ω）\n');
pr(fid, '压力系数 a_p=(1/c_f)dc_f/dp：绝热 (B/2A)/K_s = %.3e /MPa；等温（Beyer 修正）= %.3e /MPa\n', ...
    F.a_p_s*MPa, F.a_p*MPa);
% (1/c_K) dc_K/dp = [a_p − Ω/(2K_T)]/(1+Ω) − Ω/(2(1+Ω)) (dεθ/dp − dεr/dp)
% 以下全部按"每 MPa"计
term_ap   = (F.a_p*MPa)/(1 + Om);
term_rho  = -(Om/(2*F.KT/MPa))/(1 + Om);
% 几何：σθ = 5.5 p, σz = 2.52 p(闭端), σr = −0.5 p
sig_p = diag([G.sig_z_closed_per_p, G.sig_theta_per_p, G.sig_r_per_p])*MPa;   % 每 MPa 内压
eps_p = hooke(sig_p);                    % 每 MPa 内压的应变 (z, θ, r)
term_geo  = -(Om/(2*(1 + Om)))*(eps_p(2,2) - eps_p(3,3));
kK_p = term_ap + term_rho + term_geo;    % 每 MPa 内压
pr(fid, '(1/c_K)dc_K/dp [1/MPa]：a_p 项 %.3e，ρ_f 项 %.3e，几何项 %.3e，合计 %.3e\n', ...
    term_ap, term_rho, term_geo, kK_p);
kappa0 = (term_ap + term_rho)/(F.a_p*MPa);
pr(fid, 'κ(0) = 液模态压力系数 / 自由液体 a_p = %.3f（=(B/A−Ω K_s/K_T)/((1+Ω)B/A) 绝热近似 %.3f）\n', ...
    kappa0, (F.BA - Om*F.Ks/F.KT)/((1 + Om)*F.BA));
% 液模态对 σz：几何 εθ−εr = 0（各向同性泊松），只剩壁刚度声弹变化（量级 K_θ·Ω/(2(1+Ω))）
kK_sz = (Om/(2*(1 + Om)))*(2*K_nat(2));   % 环向刚度对轴向应力的声弹变化，粗估 2K_θ
pr(fid, '液模态对 σz（真实速度，物理项）≈ %.1e /MPa（壁环向刚度声弹变化粗估，几何项恒为 0）\n', kK_sz);
t0_f = G.W_E/c_K;
dt_dp_f = -t0_f*kK_p;   dt_dsz_f = -t0_f*kK_sz;
pr(fid, '延时 t0 = %.2f µs：∂t/∂p = %+.1f ns/MPa，∂t/∂σz ≈ %+.2f ns/MPa\n', t0_f*1e6, dt_dp_f*1e9, dt_dsz_f*1e9);
pr(fid, '两种液体对比：油 vs 水 的 a_p 之比（理论预测液模态压力灵敏度比）：\n');
Pw = params_pipe('7075', 'water');
Omw = 2*Pw.fluid.Ks*eps_hat;
kw = (Pw.fluid.a_p*MPa)/(1 + Omw) - (Omw/(2*Pw.fluid.KT/MPa))/(1 + Omw);
pr(fid, '  油 a_p=%.3e，水 a_p=%.3e /MPa；含管壁修正后 油 %.3e，水 %.3e /MPa，比 %.2f\n', ...
    F.a_p*MPa, Pw.fluid.a_p*MPa, (term_ap + term_rho), kw, (term_ap + term_rho)/kw);
pr(fid, '  水：Ω=%.4f  c_K=%.1f m/s  t0=%.2f µs  ∂t/∂p=%+.1f ns/MPa\n\n', ...
    Omw, Pw.fluid.c/sqrt(1 + Omw), G.W_E/(Pw.fluid.c/sqrt(1 + Omw))*1e6, ...
    -G.W_E/(Pw.fluid.c/sqrt(1 + Omw))*kw*1e9);

%% 4. 温度系数
pr(fid, '== 4. 温度系数（探头空间固定，路径长度不变）==\n');
% 壁 S0：(1/c)dc/dT = ½[(1/E)dE/dT + 2ν/(1−ν^2) dν/dT + 3α]
kT_w = 0.5*(W.dlnE_dT + 2*W.nu/(1 - W.nu^2)*W.dnu_dT + 3*W.alpha);
% 液 Korteweg：(1/c_K)dc_K/dT = a_T/(1+Ω) + Ω/(2(1+Ω)) [β_f + (1/E)dE/dT]
kT_f = F.a_T/(1 + Om) + (Om/(2*(1 + Om)))*(F.beta + W.dlnE_dT);
pr(fid, '壁：(1/c)dc/dT = %.3e /K → ∂t/∂T = %+.2f ns/K → 假应力 %.2f MPa/K\n', ...
    kT_w, -t0_w*kT_w*1e9, -t0_w*kT_w/dt_dsz);
pr(fid, '液：(1/c_K)dc_K/dT = %.3e /K → ∂t/∂T = %+.1f ns/K → 假压力 %.2f MPa/K\n', ...
    kT_f, -t0_f*kT_f*1e9, -t0_f*kT_f/dt_dp_f);
pr(fid, '（(1/E)dE/dT=%.1e、dc_f/dT=%g m/s/K 为典型值，第三列以实测为准）\n\n', W.dlnE_dT, F.dc_dT);

%% 5. 约定三：几何更新的贡献
pr(fid, '== 5. 几何更新（半径、壁厚随载荷变化）的贡献 ==\n');
fd = 0.625e3;                      % Hz·m (625 kHz × 1 mm)
cS0 = lamb_S0_phase(fd, W.cL, W.cT);
dfd = 0.01*fd;
slope = (log(lamb_S0_phase(fd + dfd, W.cL, W.cT)) - log(lamb_S0_phase(fd - dfd, W.cL, W.cT)))/(log(fd + dfd) - log(fd - dfd));
pr(fid, 'S0 在 fd=%.3f MHz·mm：c=%.1f m/s（薄板极限 %.1f），d ln c/d ln(fd) = %.4f\n', fd/1e3, cS0, c_S0, slope);
eps_z1 = hooke(diag([1 0 0])*MPa);          % 每 MPa σz 的应变
geo_w_sz = slope*eps_z1(3,3);               % 壁厚变化 → fd 变化
geo_w_p  = slope*eps_p(3,3);
pr(fid, '壁模态几何项：每 MPa σz %.2e（占 K_z %.2f%%）；每 MPa p %.2e（占 K_p %.2f%%）\n', ...
    geo_w_sz, geo_w_sz/K_true(1)*100, geo_w_p, geo_w_p/Kp_true*100);
pr(fid, '液模态几何项：每 MPa p %.2e（占 %.2f%%）；每 MPa σz 0（εθ = εr）\n', ...
    term_geo, term_geo/kK_p*100);
pr(fid, '高频（刚壁液柱）极限下半径项 (c_p^2/c_f^2 − 1)·εθ：c_p/c_f=1.1 时每 MPa p %.2e（占 %.1f%%）\n', ...
    0.21*eps_p(2,2), 0.21*eps_p(2,2)/(F.a_p*MPa)*100);
pr(fid, '结论：两族的几何项都在 1%% 以下，FEM 不必更新几何；但参考构型→真实速度的 εz 换算必须做（壁模态 %.0f%%）。\n\n', ...
    100/(W.E/MPa)/abs(K_nat(1)));

%% 6. 汇总
pr(fid, '== 6. 汇总：延时灵敏度（真实速度，标距 %g mm）==\n', G.W_E*1e3);
pr(fid, '%-22s %12s %12s %12s %10s\n', '模态', '∂t/∂σz ns/MPa', '∂t/∂p ns/MPa', '∂t/∂T ns/K', 't0 µs');
pr(fid, '%-22s %12.3f %12.3f %12.2f %10.2f\n', '壁 S0(薄板极限)', dt_dsz*1e9, dt_dp_w*1e9, -t0_w*kT_w*1e9, t0_w*1e6);
pr(fid, '%-22s %12.3f %12.1f %12.1f %10.2f\n', '液 Korteweg(长波极限)', dt_dsz_f*1e9, dt_dp_f*1e9, -t0_f*kT_f*1e9, t0_f*1e6);
S = [dt_dsz, dt_dp_w; dt_dsz_f, dt_dp_f]*1e9;
pr(fid, '2×2 矩阵 [ns/MPa]：[%8.3f %8.3f; %8.3f %8.3f]，条件数（列按各自最大值归一后）%.2f\n', ...
    S(1,1), S(1,2), S(2,1), S(2,2), cond(S./max(abs(S), [], 1)));
pr(fid, '5 ns 计时误差 → p 误差 %.3f MPa，σz 误差 %.2f MPa（其中来自 p 的传递 %.3f MPa）\n', ...
    5/abs(S(2,2)), 5/abs(S(1,1)), abs(S(1,2))/abs(S(2,2))*5/abs(S(1,1)));
fclose(fid);
fprintf('\n结果已写入 %s\n', fullfile(here, 'T1_结果.txt'));

%% 本地函数
function pr(fid, varargin)
fprintf(1, varargin{:});
fprintf(fid, varargin{:});
end
