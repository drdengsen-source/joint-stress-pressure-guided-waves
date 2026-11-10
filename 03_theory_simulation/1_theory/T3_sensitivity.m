%% T3_sensitivity.m  数量级表与图
% 用 SAFE（safe_pipe_setup / safe_pipe_solve / safe_perturb）算无预应力充液管的轴对称模态，
% 按 T2 §2.1 一阶摄动求每个模态对 单位 σz、单位 p、单位 ΔT 的波数变化，
% 换成真实速度、固定标距 W_E 的延时灵敏度（相速度版 = 频散补偿相位计时；群速度版 = 包络计时），画：
%   图1 频散（相/群速度，按液体能量占比着色）
%   图2 壁传族：∂t/∂σz、∂t/∂p(交叉：壁声弹 + 液体加载)、假应力比 随频率（上排相位版，下排群版）
%   图3 液传族：∂t/∂p、κ_p/κ_g、∂t/∂σz 随频率
%   图4 温度：两族 ∂t/∂T 随频率
%   图5 2×2 矩阵条件数与 σz 误差放大随 (f_w, f_l)（相位版）
%   图6 忽略内压的假应力偏差（相位版）
% 约定：拉为正；真实速度；单位扰动 1 MPa / 1 K；ns 计。图窗可见，PNG+FIG 存 figures\。
% 运行前提：params_pipe.m 的参数；不要在 MCP 会话里 clear all。
clear; clc;
here = fileparts(mfilename('fullpath'));
figdir = fullfile(here, 'figures');  if ~exist(figdir, 'dir'), mkdir(figdir); end
P = params_pipe('7075', 'oil');  W = P.wall;  F = P.fluid;  G = P.geom;
MPa = 1e6;
hooke = @(sig) ((1 + W.nu)*sig - W.nu*trace(sig)*eye(3))/W.E;
S = safe_pipe_setup(P, 8, 40);
D0 = tensor_to_D5(murnaghan_A(W.lam, W.mu, W.l, W.m, W.n, zeros(3), zeros(3)));
fid = fopen(fullfile(here, 'T3_结果.txt'), 'w', 'n', 'UTF-8');
pr(fid, 'T3 灵敏度随频率  %s\n管壁 %s；液体 %s（c_f=%g m/s）；W_E = %g mm\n', datestr(now, 'yyyy-mm-dd HH:MM'), W.name, F.name, F.c, G.W_E*1e3);

%% 单位扰动（T2 §1 的 (a)(b)(c)(d)）
% (1) 单位 σz = 1 MPa 均匀；液体拉回几何项；真实速度换算 +εz
sig = diag([0 0 1])*MPa;  e1 = hooke(sig);
dD_z = tensor_to_D5(murnaghan_A(W.lam, W.mu, W.l, W.m, W.n, sig, e1)) - D0;
ez = e1(3,3);  eth = e1(2,2);  tr = 2*eth + ez;
fl_z = [-F.Ks*tr, F.rho*(2*eth - tr), F.rho*(2*ez - tr)];      % [ΔK_s, Δρ_rr, Δρ_zz]
conv_z = ez;
% (2) 单位 p = 1 MPa：壁 Lamé（σz = 0），液体物理 + 拉回；换算 +εz(p)
a = G.r_i^2/(G.r_o^2 - G.r_i^2);
nG = numel(S.s.r);  dD_p = zeros(5, 5, nG);
for i = 1:nG
    r = S.s.r(i);
    sg = diag([MPa*a*(1 - G.r_o^2/r^2), MPa*a*(1 + G.r_o^2/r^2), 0]);
    dD_p(:,:,i) = tensor_to_D5(murnaghan_A(W.lam, W.mu, W.l, W.m, W.n, sg, hooke(sg))) - D0;
end
ei = hooke(diag([-MPa, MPa*a*(1 + G.r_o^2/G.r_i^2), 0]));
eth_p = ei(2,2);  ez_p = ei(3,3);  tr_p = 2*eth_p + ez_p;
dKs_ph = F.Ks*(1/F.KT + 2*F.a_p)*MPa;  drho_ph = F.rho*MPa/F.KT;      % 绝对增量
fl_p = [dKs_ph - F.Ks*tr_p, drho_ph + F.rho*(2*eth_p - tr_p), drho_ph + F.rho*(2*ez_p - tr_p)];
conv_p = ez_p;
% (3) 单位 ΔT = 1 K：物理量直接变，路径固定，不换算
dD_T = D0*W.dlnE_dT;  drs_T = -3*W.alpha*W.rho;
fl_T = [F.Ks*(2*F.a_T - F.beta), -F.rho*F.beta, -F.rho*F.beta];
% 第 4 列：单位 p 只含壁 ΔA（声弹交叉项），第 2 列减第 4 列 = 液体加载项（K_s、ρ_f 变化 + 拉回）
pert = @(md, om) [ ...
    safe_perturb(S, md, om, dD_z, fl_z(1), fl_z(2), fl_z(3), 0,     F.Ks, F.rho), ...
    safe_perturb(S, md, om, dD_p, fl_p(1), fl_p(2), fl_p(3), 0,     F.Ks, F.rho), ...
    safe_perturb(S, md, om, dD_T, fl_T(1), fl_T(2), fl_T(3), drs_T, F.Ks, F.rho), ...
    safe_perturb(S, md, om, dD_p, 0, 0, 0, 0, F.Ks, F.rho)];

%% 频率扫描
fk = sort(unique([(20:10:800)*1e3, 625e3]));  nf = numel(fk);  dl = 0.005;   % 网格含 625 kHz 工作点
cmin = 900;  cmax = 30000;
R = struct('f', {}, 'cp', {}, 'cg', {}, 'eta', {}, 'k', {}, 'dk', {}, 'ddk', {});
tic;
for n = 1:nf
    om = 2*pi*fk(n);
    md  = safe_pipe_solve(S, om,          D0, W.rho, F.Ks, F.rho, cmin, cmax);
    mdp = safe_pipe_solve(S, om*(1 + dl), D0, W.rho, F.Ks, F.rho, cmin, cmax);
    mdm = safe_pipe_solve(S, om*(1 - dl), D0, W.rho, F.Ks, F.rho, cmin, cmax);
    nm = numel(md);
    rec = struct('f', fk(n), 'cp', [md.cp], 'cg', [md.cg], 'eta', [md.eta_f], 'k', [md.k], ...
                 'dk', zeros(nm, 4), 'ddk', nan(nm, 4));
    for m = 1:nm
        rec.dk(m,:) = pert(md(m), om);
        % 群速度版：ω(1±δ) 处按最近相速度配对，相对差 >3% 视为跨模态，弃
        [dp_, jp] = min(abs([mdp.cp] - md(m).cp)/md(m).cp);
        [dm_, jm] = min(abs([mdm.cp] - md(m).cp)/md(m).cp);
        if ~isempty(jp) && ~isempty(jm) && dp_ < 0.03 && dm_ < 0.03
            rec.ddk(m,:) = (pert(mdp(jp), om*(1 + dl)) - pert(mdm(jm), om*(1 - dl)))/(2*dl*om);
        end
    end
    R(n) = rec;
    if mod(n, 10) == 0, fprintf('  %d/%d  %.0f s\n', n, nf, toc); end
end

%% 换成灵敏度
% 相速度版：K^true = −Δk/k + conv；S = −t0p·K^true（ns/MPa 或 ns/K）
% 群速度版：Δt = W_E·dΔk/dω − t0g·conv
conv = [conv_z, conv_p, 0, 0];                  % 第 4 列（壁 ΔA 分项）只报自然量
allF = []; allCP = []; allCG = []; allETA = []; allSp = []; allSg = []; allK = [];
for n = 1:nf
    rec = R(n);  nm = numel(rec.cp);
    Knat = -rec.dk./rec.k(:);                   % nm×4，每 MPa / 每 K
    Ktrue = Knat + conv;                        % 液模态含拉回 −εz，换算后抵消
    t0p = G.W_E./rec.cp(:);  t0g = G.W_E./rec.cg(:);
    Sp = -t0p.*Ktrue*1e9;                       % ns
    Sg = (G.W_E*rec.ddk - t0g.*conv)*1e9;
    allF = [allF; rec.f*ones(nm,1)]; allCP = [allCP; rec.cp(:)]; allCG = [allCG; rec.cg(:)];
    allETA = [allETA; rec.eta(:)]; allSp = [allSp; Sp]; allSg = [allSg; Sg]; allK = [allK; Ktrue];
end
% 族划分
isWall = allCP > 4300 & allCP < 6500 & allETA < 0.5;
isLiq  = allCP > 1000 & allCP < 2600 & allETA > 0.5;
% 壁传族每频率取壁能量占比最大者（无 η_f<0.5 的模态 = 耦合区，留空）
wallIdx = nan(nf, 1);
for n = 1:nf
    ii = find(allF == fk(n) & isWall);
    if ~isempty(ii), [~, j] = min(allETA(ii)); wallIdx(n) = ii(j); end
end
wallIdx = wallIdx(~isnan(wallIdx));  fw = allF(wallIdx);  etaW = allETA(wallIdx);
% 干净窗口：壁模态 η_f < 0.10 的连续频段
clean = etaW(:) < 0.10;  edges = diff([0; clean; 0]);
wstart = fw(edges(1:end-1) == 1);  wend = fw(edges(2:end) == -1);

%% T1 对照值
A00 = murnaghan_A(W.lam, W.mu, W.l, W.m, W.n, zeros(3), zeros(3));
KzT1 = (s0_reduced_modulus(murnaghan_A(W.lam, W.mu, W.l, W.m, W.n, diag([1 0 0])*MPa, hooke(diag([1 0 0])*MPa))) - s0_reduced_modulus(A00))/(2*s0_reduced_modulus(A00)) + 1/(W.E/MPa);
KthT1 = (s0_reduced_modulus(murnaghan_A(W.lam, W.mu, W.l, W.m, W.n, diag([0 1 0])*MPa, hooke(diag([0 1 0])*MPa))) - s0_reduced_modulus(A00))/(2*s0_reduced_modulus(A00)) - W.nu/(W.E/MPa);
KrT1  = (s0_reduced_modulus(murnaghan_A(W.lam, W.mu, W.l, W.m, W.n, diag([0 0 1])*MPa, hooke(diag([0 0 1])*MPa))) - s0_reduced_modulus(A00))/(2*s0_reduced_modulus(A00)) - W.nu/(W.E/MPa);
KpT1 = G.sig_theta_per_p*KthT1 + G.sig_r_per_p*KrT1;
t0T1 = G.W_E/W.c_plate;
S11T1 = -t0T1*KzT1*1e9;  S12T1 = -t0T1*KpT1*1e9;
eps_hat = ((G.r_o^2 + G.r_i^2)/(G.r_o^2 - G.r_i^2) + W.nu)/W.E;  Om = 2*F.Ks*eps_hat;
cK = F.c/sqrt(1 + Om);  kK = ((F.a_p*MPa) - Om/(2*F.KT/MPa))/(1 + Om);
S22T1 = -(G.W_E/cK)*kK*1e9;
ap = F.a_p*MPa;

%% 图 1 频散
fig1 = figure('Visible', 'on', 'Name', 'T3-1 dispersion', 'Position', [60 60 1100 450], 'Color', 'w');
for pnl = 1:2
    subplot(1, 2, pnl); hold on; box on;
    if pnl == 1, y = allCP; else, y = allCG; end
    scatter(allF/1e3, y, 12, allETA, 'filled');
    colormap(turbo); cb = colorbar; cb.Label.String = 'fluid energy fraction \eta_f';
    yline(F.c, '--k', 'c_f'); yline(cK, ':k', 'c_K'); yline(W.c_plate, '--k', 'c_{plate}'); yline(W.c_bar, ':k', 'c_{bar}');
    xline(625, '-.r'); xline(250, '-.b'); xline(400, '-.b');
    xlabel('f (kHz)');
    if pnl == 1, ylabel('phase velocity (m/s)'); title('(a) c_p(f)'); else, ylabel('group velocity (m/s)'); title('(b) c_g(f)'); end
    ylim([800 7000]); xlim([0 800]);
end
exportgraphics(fig1, fullfile(figdir, 'T3_1_dispersion.png'), 'Resolution', 200);
savefig(fig1, fullfile(figdir, 'T3_1_dispersion.fig'));

%% 图 2 壁传族（上排相位版，下排群版）
fig2 = figure('Visible', 'on', 'Name', 'T3-2 wall family', 'Position', [80 80 1250 720], 'Color', 'w');
subplot(2,3,1); hold on; box on; shade(wstart, wend);
plot(fw/1e3, allSp(wallIdx,1), 'b-', 'LineWidth', 1.3); yline(S11T1, '--k', 'T1 thin plate');
xlabel('f_w (kHz)'); ylabel('\partial t/\partial\sigma_z (ns/MPa)'); title('(a) phase: axial stress'); xlim([0 800]); ylim([0.8 1.4]);
subplot(2,3,2); hold on; box on; shade(wstart, wend);
h1 = plot(fw/1e3, allSp(wallIdx,2), 'b-', 'LineWidth', 1.3);
h2 = plot(fw/1e3, allSp(wallIdx,4), 'r:', 'LineWidth', 1.1);
h3 = plot(fw/1e3, allSp(wallIdx,2) - allSp(wallIdx,4), 'r--', 'LineWidth', 1.1);
yline(S12T1, '--k', 'T1 thin plate');
xlabel('f_w (kHz)'); ylabel('\partial t/\partial p (ns/MPa)'); title('(b) phase: pressure cross-term');
legend([h1 h2 h3], {'total', 'wall acoustoelastic', 'fluid loading'}, 'Location', 'southwest'); xlim([0 800]); ylim([-35 5]);
subplot(2,3,3); hold on; box on; shade(wstart, wend);
plot(fw/1e3, allSp(wallIdx,2)./allSp(wallIdx,1), 'b-', 'LineWidth', 1.3); yline(KpT1/KzT1, '--k', 'T1');
xlabel('f_w (kHz)'); ylabel('apparent \sigma_z per MPa of p'); title('(c) phase: bias of single-mode method'); xlim([0 800]); ylim([-30 2]);
subplot(2,3,4); hold on; box on; shade(wstart, wend);
plot(fw/1e3, allSg(wallIdx,1), 'r-', 'LineWidth', 1.3); yline(S11T1, '--k');
xlabel('f_w (kHz)'); ylabel('\partial t/\partial\sigma_z (ns/MPa)'); title('(d) group: axial stress (clipped)'); xlim([0 800]); ylim([-5 8]);
subplot(2,3,5); hold on; box on; shade(wstart, wend);
plot(fw/1e3, allSg(wallIdx,2), 'r-', 'LineWidth', 1.3); yline(S12T1, '--k');
xlabel('f_w (kHz)'); ylabel('\partial t/\partial p (ns/MPa)'); title('(e) group: pressure cross-term (clipped)'); xlim([0 800]); ylim([-300 300]);
subplot(2,3,6); hold on; box on; shade(wstart, wend);
plot(fw/1e3, allSg(wallIdx,2)./allSg(wallIdx,1), 'r-', 'LineWidth', 1.3); yline(KpT1/KzT1, '--k');
xlabel('f_w (kHz)'); ylabel('apparent \sigma_z per MPa of p'); title('(f) group: bias (clipped)'); xlim([0 800]); ylim([-300 300]);
exportgraphics(fig2, fullfile(figdir, 'T3_2_wall_family.png'), 'Resolution', 200);
savefig(fig2, fullfile(figdir, 'T3_2_wall_family.fig'));

%% 图 3 液传族
fig3 = figure('Visible', 'on', 'Name', 'T3-3 liquid family', 'Position', [100 100 1250 420], 'Color', 'w');
li = find(isLiq);
kap_p = allK(li,2)/ap;
kap_g = (-allCG(li).*(allSg(li,2)*1e-9 + (G.W_E./allCG(li))*conv_p)/G.W_E + conv_p)/ap;
subplot(1,3,1); hold on; box on;
scatter(allF(li)/1e3, allSp(li,2), 14, allETA(li), 'filled'); colormap(turbo); caxis([0.5 1]);
scatter(allF(li)/1e3, allSg(li,2), 10, 'k');
yline(S22T1, '--k', 'Korteweg (thick wall)');
xlabel('f_l (kHz)'); ylabel('\partial t/\partial p (ns/MPa)'); title('(a) liquid modes: pressure (colour = phase, black = group)'); xlim([0 800]); ylim([-800 400]);
subplot(1,3,2); hold on; box on;
scatter(allF(li)/1e3, kap_p, 14, 'b', 'filled'); scatter(allF(li)/1e3, kap_g, 14, 'r', 'filled');
yline(kK/ap, '--k', '\kappa(0)');
xlabel('f_l (kHz)'); ylabel('\kappa'); title('(b) \kappa_p (blue), \kappa_g (red)'); xlim([0 800]); ylim([0 2]);
subplot(1,3,3); hold on; box on;
scatter(allF(li)/1e3, allSp(li,1), 14, allETA(li), 'filled'); caxis([0.5 1]);
scatter(allF(li)/1e3, allSg(li,1), 10, 'k');
xlabel('f_l (kHz)'); ylabel('\partial t/\partial\sigma_z (ns/MPa)'); title('(c) liquid modes: axial-stress cross-talk'); xlim([0 800]); ylim([-8 8]);
exportgraphics(fig3, fullfile(figdir, 'T3_3_liquid_family.png'), 'Resolution', 200);
savefig(fig3, fullfile(figdir, 'T3_3_liquid_family.fig'));

%% 图 4 温度
fig4 = figure('Visible', 'on', 'Name', 'T3-4 temperature', 'Position', [120 120 900 400], 'Color', 'w');
subplot(1,2,1); hold on; box on; shade(wstart, wend);
h4 = plot(fw/1e3, allSp(wallIdx,3), 'b-', 'LineWidth', 1.2);
h5 = plot(fw/1e3, allSg(wallIdx,3), 'r-', 'LineWidth', 1.2);
xlabel('f_w (kHz)'); ylabel('\partial t/\partial T (ns/K)'); title('(a) wall mode'); legend([h4 h5], {'phase', 'group'}, 'Location', 'best'); xlim([0 800]); ylim([-50 60]);
subplot(1,2,2); hold on; box on;
scatter(allF(li)/1e3, allSp(li,3), 14, allETA(li), 'filled'); colormap(turbo); caxis([0.5 1]);
scatter(allF(li)/1e3, allSg(li,3), 10, 'k');
xlabel('f_l (kHz)'); ylabel('\partial t/\partial T (ns/K)'); title('(b) liquid modes (colour = phase, black = group)'); xlim([0 800]); ylim([-300 600]);
exportgraphics(fig4, fullfile(figdir, 'T3_4_temperature.png'), 'Resolution', 200);
savefig(fig4, fullfile(figdir, 'T3_4_temperature.fig'));

%% 液传族代表模态：每个 f_l 取 η_f>0.6、c_g∈[1250,1500]（~70–80 µs 到达的包）中 η_f 最大者
fl_grid = fk(fk >= 150e3 & fk <= 500e3);
Lsel = nan(size(fl_grid));
for n = 1:numel(fl_grid)
    ii = find(allF == fl_grid(n) & isLiq & allETA > 0.6 & allCG > 1250 & allCG < 1500);
    if ~isempty(ii), [~, j] = max(allETA(ii)); Lsel(n) = ii(j); end
end

%% 图 5 条件数与误差放大（相位版）
fw_grid = fw(fw >= 300e3);
condMap = nan(numel(fl_grid), numel(fw_grid)); errMap = condMap;
for i = 1:numel(fl_grid)
    if isnan(Lsel(i)), continue; end
    for j = 1:numel(fw_grid)
        w_ = wallIdx(fw == fw_grid(j));
        Sm = [allSp(w_,1), allSp(w_,2); allSp(Lsel(i),1), allSp(Lsel(i),2)];
        condMap(i,j) = cond(Sm./max(abs(Sm), [], 1));
        errMap(i,j) = sqrt(Sm(2,2)^2 + Sm(1,2)^2)/abs(det(Sm));   % σ_σz 每 1 ns 计时误差 (MPa/ns)
    end
end
fig5 = figure('Visible', 'on', 'Name', 'T3-5 condition & error', 'Position', [140 140 1000 420], 'Color', 'w');
subplot(1,2,1); imagesc(fw_grid/1e3, fl_grid/1e3, condMap, [1 3]); axis xy; colorbar; xlabel('f_w (kHz)'); ylabel('f_l (kHz)'); title('(a) cond(S), columns normalised (phase)');
subplot(1,2,2); imagesc(fw_grid/1e3, fl_grid/1e3, errMap, [0.7 1.3]); axis xy; colorbar; xlabel('f_w (kHz)'); ylabel('f_l (kHz)'); title('(b) \sigma_{\sigma z} per 1 ns timing error (MPa)');
colormap(fig5, parula);
exportgraphics(fig5, fullfile(figdir, 'T3_5_condition_error.png'), 'Resolution', 200);
savefig(fig5, fullfile(figdir, 'T3_5_condition_error.fig'));

%% 推荐工作点：350–800 kHz 内壁模态 η_f 最小处；625 kHz 作对照
inr = fw >= 350e3;  [~, jj] = min(etaW(inr));  tmp = find(inr);  jbest = tmp(jj);  fw_best = fw(jbest);
[~, j625] = min(abs(fw - 625e3));
[~, j500] = min(abs(fw - 500e3));
biasP = allSp(wallIdx,2)./allSp(wallIdx,1);

%% 图 6 忽略内压的偏差（相位版）
pv = 0:0.5:20;
fig6 = figure('Visible', 'on', 'Name', 'T3-6 pressure bias', 'Position', [160 160 1000 420], 'Color', 'w');
subplot(1,2,1); hold on; box on; shade(wstart, wend);
plot(fw/1e3, biasP, 'b-', 'LineWidth', 1.3); yline(KpT1/KzT1, '--k', 'T1 thin plate');
xlabel('f_w (kHz)'); ylabel('apparent \sigma_z per MPa of p (MPa/MPa)'); title('(a) bias coefficient vs f_w (phase timing)'); xlim([0 800]); ylim([-30 2]);
subplot(1,2,2); hold on; box on;
plot(pv, biasP(jbest)*pv, 'b-', pv, biasP(j625)*pv, 'r-', pv, (KpT1/KzT1)*pv, '--k', 'LineWidth', 1.3);
yline(-3.8, ':k', 'old claim \pm3.8 MPa'); yline(3.8, ':k');
xlabel('internal pressure p (MPa)'); ylabel('apparent axial stress error (MPa)'); title('(b) bias of single-mode method');
legend(sprintf('%.0f kHz (recommended)', fw_best/1e3), '625 kHz', 'T1 thin plate', 'Location', 'southwest');
exportgraphics(fig6, fullfile(figdir, 'T3_6_pressure_bias.png'), 'Resolution', 200);
savefig(fig6, fullfile(figdir, 'T3_6_pressure_bias.fig'));

%% 数值表
pr(fid, '\n== 干净窗口（壁模态 η_f < 0.10 的连续频段，c_f=%g 时；截止频率随实测 c_f 成比例移动）==\n', F.c);
for k = 1:numel(wstart), pr(fid, '  %.0f–%.0f kHz（中心 %.0f）\n', wstart(k)/1e3, wend(k)/1e3, (wstart(k)+wend(k))/2e3); end
pr(fid, '推荐壁模态工作点：%.0f kHz（350 kHz 以上 η_f 最小 = %.3f）；625 kHz 处 η_f = %.3f\n', fw_best/1e3, etaW(jbest), etaW(j625));
pr(fid, '\n== 壁传族（每频率 η_f<0.5 中壁能量占比最大的模态；真实速度；ns/MPa、ns/K）==\n');
pr(fid, '%7s %8s %8s %6s | %8s %8s %8s %8s %8s | %8s %8s %8s | %8s\n', 'f/kHz', 'cp', 'cg', 'eta_f', 'S11 ph', 'S12 ph', 'S12wall', 'S12fld', 'S13 ph', 'S11 gr', 'S12 gr', 'S13 gr', 'bias ph');
for n = 1:numel(wallIdx)
    i = wallIdx(n);
    if mod(round(fw(n)/1e3), 50) == 0 || fw(n) == 625e3
        pr(fid, '%7.0f %8.1f %8.1f %6.3f | %8.3f %8.3f %8.3f %8.3f %8.2f | %8.3f %8.1f %8.1f | %8.3f\n', fw(n)/1e3, allCP(i), allCG(i), allETA(i), ...
            allSp(i,1), allSp(i,2), allSp(i,4), allSp(i,2) - allSp(i,4), allSp(i,3), allSg(i,1), allSg(i,2), allSg(i,3), biasP(n));
    end
end
pr(fid, 'T1 薄板极限：S11 = %.3f，S12 = %.3f ns/MPa（仅声弹交叉项），bias = %.3f\n', S11T1, S12T1, KpT1/KzT1);
pr(fid, '\n== 液传族：250–400 kHz 所有 η_f>0.5 且 c_g∈[1100,1700] 的模态（相位版 / 群版）==\n');
pr(fid, '%7s %8s %8s %6s | %8s %8s %8s | %8s %8s %8s | %6s %6s\n', 'f/kHz', 'cp', 'cg', 'eta_f', 'S21 ph', 'S22 ph', 'S23 ph', 'S21 gr', 'S22 gr', 'S23 gr', 'kap_p', 'kap_g');
for f0 = [250e3 300e3 350e3 400e3]
    ii = find(allF == f0 & isLiq & allCG > 1100 & allCG < 1700);
    for i = ii'
        kg = (-allCG(i)*(allSg(i,2)*1e-9 + (G.W_E/allCG(i))*conv_p)/G.W_E + conv_p)/ap;
        pr(fid, '%7.0f %8.1f %8.1f %6.3f | %8.3f %8.1f %8.1f | %8.3f %8.1f %8.1f | %6.3f %6.3f\n', f0/1e3, allCP(i), allCG(i), allETA(i), ...
            allSp(i,1), allSp(i,2), allSp(i,3), allSg(i,1), allSg(i,2), allSg(i,3), allK(i,2)/ap, kg);
    end
end
pr(fid, 'T1 厚壁 Korteweg：S22 = %.1f ns/MPa，κ(0) = %.3f\n', S22T1, kK/ap);
% 工作点矩阵（相位版）
[~, n300] = min(abs(fl_grid - 300e3));
for jw = [jbest, j500, j625]
    if isnan(Lsel(n300)), break; end
    iw = wallIdx(jw);  il = Lsel(n300);
    Sm = [allSp(iw,1), allSp(iw,2); allSp(il,1), allSp(il,2)];
    pr(fid, '\n工作点 (%.0f kHz 壁, 300 kHz 液；相位计时)：S = [%.3f %.3f; %.3f %.1f] ns/MPa，cond(列归一) = %.2f\n', fw(jw)/1e3, Sm(1,1), Sm(1,2), Sm(2,1), Sm(2,2), cond(Sm./max(abs(Sm), [], 1)));
    pr(fid, '  1 ns 计时误差 → σ_p = %.4f MPa，σ_σz = %.3f MPa；5 ns → σ_σz = %.2f MPa\n', sqrt(Sm(1,1)^2 + Sm(2,1)^2)/abs(det(Sm)), sqrt(Sm(2,2)^2 + Sm(1,2)^2)/abs(det(Sm)), 5*sqrt(Sm(2,2)^2 + Sm(1,2)^2)/abs(det(Sm)));
    pr(fid, '  温度：0.1 K 误差 → σz %.2f MPa（S13=%.2f ns/K），p %.3f MPa（S23=%.1f ns/K）\n', 0.1*abs(allSp(iw,3)/Sm(1,1)), allSp(iw,3), 0.1*abs(allSp(il,3)/Sm(2,2)), allSp(il,3));
    pr(fid, '  忽略内压（单模态法）：20 MPa 内压 → 假应力 %.1f MPa（T1 薄板 %.1f）\n', 20*biasP(jw), 20*KpT1/KzT1);
end
fclose(fid);
save(fullfile(here, 'T3_sweep.mat'), 'R', 'fk', 'allF', 'allCP', 'allCG', 'allETA', 'allSp', 'allSg', 'allK', 'wallIdx', 'fw', 'conv', 'P', 'wstart', 'wend', 'fw_best');
fprintf('\n完成：结果 T3_结果.txt，数据 T3_sweep.mat，图 %s\n', figdir);
figure(fig1); figure(fig2); figure(fig3); figure(fig4); figure(fig5); figure(fig6); drawnow;

%% 本地函数
function pr(fid, varargin)
fprintf(1, varargin{:});
fprintf(fid, varargin{:});
end

function shade(ws, we)
% 干净窗口淡绿底
for k = 1:numel(ws)
    patch([ws(k) we(k) we(k) ws(k)]/1e3, [-1e9 -1e9 1e9 1e9], [0.9 1 0.9], 'EdgeColor', 'none');
end
end
