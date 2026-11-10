%% T3_check_perturb.m  摄动公式的数值校核（T2 §8 第 3、4 条）
%   1. 摄动 Δk vs 直接法（把 ΔA、ΔK_s、Δρ_f 放进 SAFE 重解）——验证一阶摄动与线性度
%   2. 低频慢波对 p 的灵敏度 → T1 厚壁 Korteweg 闭式
%   3. 高频快波对 σz 的灵敏度 → T1 薄板 S0 闭式（差异 = 有限 fd + 曲率 + 液体耦合）
% 运行：MATLAB 里 cd 到 1_theory 后 T3_check_perturb
clear; clc;
P = params_pipe('7075','oil'); W = P.wall; F = P.fluid; G = P.geom;
MPa = 1e6;  hooke = @(sig) ((1+W.nu)*sig - W.nu*trace(sig)*eye(3))/W.E;
S = safe_pipe_setup(P, 8, 40);
D0 = tensor_to_D5(murnaghan_A(W.lam, W.mu, W.l, W.m, W.n, zeros(3), zeros(3)));
cmin = 900; cmax = 30000;

%% 1. 单位 σz：摄动 vs 直接法（只动壁，液体拉回项两边都不加）
sig = diag([0 0 1])*MPa;  eps1 = hooke(sig);
dD_z = tensor_to_D5(murnaghan_A(W.lam, W.mu, W.l, W.m, W.n, sig, eps1)) - D0;
fprintf('== 1. 单位 σz：摄动 vs 直接法（壁 ΔA，100 MPa 直接解）==\n');
fprintf('%8s %9s %7s %13s %13s %9s\n', 'f/kHz', 'cp', 'eta_f', 'Knat 摄动', 'Knat 直接', '差%');
for f = [20e3 300e3 625e3 800e3]
    om = 2*pi*f;
    md  = safe_pipe_solve(S, om, D0,          W.rho, F.Ks, F.rho, cmin, cmax);
    md2 = safe_pipe_solve(S, om, D0 + 100*dD_z, W.rho, F.Ks, F.rho, cmin, cmax);
    for m = 1:numel(md)
        dk = safe_perturb(S, md(m), om, dD_z, 0, 0, 0, 0, F.Ks, F.rho);
        Kpert = -dk/md(m).k;
        [dc, j] = min(abs([md2.cp] - md(m).cp)/md(m).cp);
        Kdir = (md(m).k/md2(j).k - 1)/100;
        if dc < 0.05
            fprintf('%8.0f %9.1f %7.3f %13.4e %13.4e %9.2f\n', f/1e3, md(m).cp, md(m).eta_f, Kpert, Kdir, (Kdir/Kpert - 1)*100);
        end
    end
end

%% 2. 单位 p：液体物理项（ΔK_s、Δρ_f 各向同性）+ 壁 Lamé ΔA，摄动 vs 直接法（20 MPa）
a = G.r_i^2/(G.r_o^2 - G.r_i^2);
nG = numel(S.s.r);  dD_p = zeros(5,5,nG);
for i = 1:nG
    r = S.s.r(i);
    sr = MPa*a*(1 - G.r_o^2/r^2);  sth = MPa*a*(1 + G.r_o^2/r^2);
    sg = diag([sr sth 0]);
    dD_p(:,:,i) = tensor_to_D5(murnaghan_A(W.lam, W.mu, W.l, W.m, W.n, sg, hooke(sg))) - D0;
end
dKs = F.Ks*(1/F.KT + 2*F.a_p)*MPa;  drho = F.rho*MPa/F.KT;     % 绝对密度增量 Δρ_f = ρ_f Δp/K_T
stp = 2;                                  % 直接法步长 MPa（20 MPa 会跨模态配对失败）
fprintf('\n== 2. 单位 p：摄动 vs 直接法（壁 Lamé + 液体物理项，%g MPa 直接解）==\n', stp);
fprintf('%8s %9s %7s %13s %13s %9s | %11s %11s %11s\n', 'f/kHz', 'cp', 'eta_f', 'Knat 摄动', 'Knat 直接', '差%', '壁ΔA项', '液K_s项', '液ρ_f项');
for f = [5e3 20e3 300e3 625e3]
    om = 2*pi*f;
    md  = safe_pipe_solve(S, om, D0,    W.rho, F.Ks,         F.rho,         cmin, cmax);
    D_ps = D0 + stp*dD_p;
    md2 = safe_pipe_solve(S, om, D_ps, W.rho, F.Ks + stp*dKs, F.rho + stp*drho, cmin, cmax);
    for m = 1:numel(md)
        dk  = safe_perturb(S, md(m), om, dD_p, dKs, drho, drho, 0, F.Ks, F.rho);
        dkw = safe_perturb(S, md(m), om, dD_p, 0, 0, 0, 0, F.Ks, F.rho);
        dkK = safe_perturb(S, md(m), om, zeros(5), dKs, 0, 0, 0, F.Ks, F.rho);
        dkr = safe_perturb(S, md(m), om, zeros(5), 0, drho, drho, 0, F.Ks, F.rho);
        Kpert = -dk/md(m).k;
        [dc, j] = min(abs([md2.cp] - md(m).cp)/md(m).cp);
        Kdir = (md(m).k/md2(j).k - 1)/stp;
        if dc < 0.02
            fprintf('%8.0f %9.1f %7.3f %13.4e %13.4e %9.2f | %11.3e %11.3e %11.3e\n', f/1e3, md(m).cp, md(m).eta_f, Kpert, Kdir, (Kdir/Kpert - 1)*100, -dkw/md(m).k, -dkK/md(m).k, -dkr/md(m).k);
        end
    end
end

%% 3. 极限对照
eps_hat = ((G.r_o^2 + G.r_i^2)/(G.r_o^2 - G.r_i^2) + W.nu)/W.E;
Om = 2*F.Ks*eps_hat;
kK = ((F.a_p*MPa) - Om/(2*F.KT/MPa))/(1 + Om);           % 厚壁 Korteweg，不含几何项与壁刚度项
om = 2*pi*5e3;
md = safe_pipe_solve(S, om, D0, W.rho, F.Ks, F.rho, cmin, cmax);
[~, m] = min([md.cp]);
dk  = safe_perturb(S, md(m), om, zeros(5), dKs, drho, drho, 0, F.Ks, F.rho);   % 只含液体物理项，与闭式同口径
dkw = safe_perturb(S, md(m), om, dD_p, 0, 0, 0, 0, F.Ks, F.rho);               % 壁 Lamé 声弹项（闭式最后一项）
fprintf('\n== 3. 极限对照 ==\n');
fprintf('慢波 5 kHz：cp=%.1f（厚壁 Korteweg %.1f）；(1/c)dc/dp 液体项 摄动 %.4e，闭式 %.4e /MPa，差 %.2f%%；壁刚度项 %.3e\n', ...
    md(m).cp, F.c/sqrt(1+Om), -dk/md(m).k, kK, ((-dk/md(m).k)/kK - 1)*100, -dkw/md(m).k);
% 快波高 fd：与薄板 S0 闭式比（T1：K_z^nat = −7.607e−5）
A1 = murnaghan_A(W.lam, W.mu, W.l, W.m, W.n, diag([1 0 0])*MPa, hooke(diag([1 0 0])*MPa));   % 板坐标 1=z
A00 = murnaghan_A(W.lam, W.mu, W.l, W.m, W.n, zeros(3), zeros(3));
KzT1 = (s0_reduced_modulus(A1) - s0_reduced_modulus(A00))/(2*s0_reduced_modulus(A00));
for f = [625e3 800e3]
    om = 2*pi*f;
    md = safe_pipe_solve(S, om, D0, W.rho, F.Ks, F.rho, 4300, 6500);
    [~, m] = min([md.eta_f]);
    dk = safe_perturb(S, md(m), om, dD_z, 0, 0, 0, 0, F.Ks, F.rho);
    fprintf('快波 %g kHz：cp=%.1f eta_f=%.3f；K_z^nat 摄动 %.4e，T1 薄板 %.4e /MPa，差 %.2f%%\n', ...
        f/1e3, md(m).cp, md(m).eta_f, -dk/md(m).k, KzT1, ((-dk/md(m).k)/KzT1 - 1)*100);
end
% 空管（去掉液体：K_s→极小，ρ_f→极小 不可行），改用板极限：fd→大时曲率影响小；此处只报告差异，不强求一致
