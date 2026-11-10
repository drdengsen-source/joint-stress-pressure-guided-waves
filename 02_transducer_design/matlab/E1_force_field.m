function E1_force_field()
%E1_FORCE_FIELD  E1 力场方案：环形永磁体在铝管外壁处的静磁场与 EMAT 洛伦兹力分布（解析）
% 输入来源：
%   初始方案：磁体 N45H 轴向宽 8、径向厚 2、内半径 8 mm，每探头两块、中心距 8.7 mm；
%       线圈 0.15 mm 漆包线 29 匝单层（宽 4.35 mm）；目标 λ_w = 10.4 mm，液模态 λ = 6.3 / 5.2 / 4.5 mm
%   1_theory\params_pipe.m：管外半径 r_o = 6.5 mm
%   N45H 剩磁取 1.35 T（牌号范围 1.32–1.38 T）
% 方法：均匀磁化 → 等效面电流 → 圆环电流场椭圆积分解析式（loopB.m、ringB.m）。周围介质 μ_r = 1，
%   解析解就是静磁 FEM 的精确参考。涡流按良导体镜像：K_φ(z) = −2 H_z,coil(r_o, z)（coilK.m）。
%   洛伦兹力面密度 f = K × B：f_z = −K_φ B_r（轴向），f_r = +K_φ B_z（径向）；乘周长 2π r_o 得每单位长度
%   的力，单位 N/m 每安培线圈电流。空间谱 F(λ) = |∫ f(z) e^{−i 2π z/λ} dz|，单位 N 每安培。
% 符号约定：z 沿管轴，r 径向；径向充磁 sgn = +1 为向外（内孔面 S 极），轴向充磁 sgn = +1 沿 +z。
%   图中磁体标签：R+ / R− 径向向外 / 向内，A→ / A← 轴向 +z / −z；线圈黑条，+/− 为电流方向。
% 单位：SI 计算；图轴 mm、T、N/(m·A)、N/A；文本结果 E1_结果.txt。
% 输出：figures\E1_1_fields.*、E1_2_forces.*、E1_3_dual.*、E1_4_variants.*；E1_结果.txt；E1_fields.mat（给 E2/E3）
% 运行：在本目录 matlab -batch "E1_force_field"，约 1 min。

here   = fileparts(mfilename('fullpath'));
figdir = fullfile(here, 'figures');
if ~exist(figdir, 'dir'), mkdir(figdir); end
addpath(here);
addpath(fullfile(here, '..', '1_theory'));
P = params_pipe('7075', 'oil');

FONT = 'TimesSimSun';
try
    if any(strcmpi(FONT, listfonts))
        set(groot, 'defaultAxesFontName', FONT, 'defaultTextFontName', FONT, 'defaultLegendFontName', FONT);
    end
catch
end
set(groot, 'defaultTextInterpreter', 'tex', 'defaultAxesTickLabelInterpreter', 'tex', ...
    'defaultLegendInterpreter', 'tex', 'defaultAxesFontSize', 9);

%% 1. 参数
mu0  = 4e-7*pi;
Brem = 1.35;                                  % N45H 剩磁 T
geo.rsurf = P.geom.r_o;                       % 管外壁 6.5 mm（涡流层 δ ≤ 0.21 mm，取外壁面）
geo.M     = Brem/mu0;                         % 磁化强度 A/m
geo.h     = 0.5e-3;                           % 线圈导线中心离壁
geo.rc    = geo.rsurf + geo.h;
geo.pitch = 0.15e-3;                          % 0.15 mm 漆包线单层匝距
geo.z     = (-32:0.02:32)'*1e-3;              % 评估轴（列）
geo.lam   = (2:0.02:30)*1e-3;                 % 空间谱 λ 轴（行）
geo.lam_t = [10.4 6.3 5.2 4.5]*1e-3;          % 目标波长：壁模态 520 kHz；液模态 300 kHz 三候选
mgE = struct('ri', 8.0e-3, 'ro', 10.0e-3, 'W', 8.0e-3);   % 现有 N45H 环
d0 = 8.7e-3;                                  % 现有中心距
d1 = 10.4e-3;                                 % 拟改中心距 = λ_w

mk = @(type, zc, sgn, mg) struct('type', type, 'zc', zc, 'sgn', sgn, 'mg', mg);
ck = @(zc, Nt, sgn) struct('zc', zc, 'Nt', Nt, 'sgn', sgn);

%% 2. 自检：单环轴上场、轴向环轴上场闭式
zt = (-20:0.5:20)'*1e-3;  a0 = 9e-3;
[~, bz] = loopB(a0, 0, 1, 1e-9, zt);
err1 = max(abs(bz - mu0*a0^2./(2*(a0^2 + zt.^2).^1.5)))/max(abs(bz));
[~, bzA] = ringB('A', mgE, 0, +1, geo.M, 1e-9, zt);
g  = @(R) (zt + mgE.W/2)./sqrt((zt + mgE.W/2).^2 + R^2) - (zt - mgE.W/2)./sqrt((zt - mgE.W/2).^2 + R^2);
bzA_ref = Brem/2*(g(mgE.ro) - g(mgE.ri));
err2 = max(abs(bzA - bzA_ref))/max(abs(bzA_ref));
fprintf('自检：单环轴上场相对误差 %.1e；轴向环轴上场相对误差 %.1e\n', err1, err2);

%% 3. 现有磁体的七种布置
A = {};
A{end+1} = arrangement('a0', '(a) 现状 径向环++ d=8.7 磁下', ...
    [mk('R',-d0/2,+1,mgE), mk('R',+d0/2,+1,mgE)], [ck(-d0/2,29,+1), ck(+d0/2,29,+1)], geo);
A{end+1} = arrangement('a1', '(a) 径向环++ d=10.4 磁下', ...
    [mk('R',-d1/2,+1,mgE), mk('R',+d1/2,+1,mgE)], [ck(-d1/2,29,+1), ck(+d1/2,29,+1)], geo);
A{end+1} = arrangement('b1', '(b) 径向环+− d=10.4 缝内 13 匝', ...
    [mk('R',-d1/2,+1,mgE), mk('R',+d1/2,-1,mgE)], ck(0,13,+1), geo);
A{end+1} = arrangement('c1', '(c) 轴向环 环中央 29 匝', ...
    mk('A',0,+1,mgE), ck(0,29,+1), geo);
A{end+1} = arrangement('c2', '(c′) 轴向环 曲折 2.6 mm 16 匝×3', ...
    mk('A',0,+1,mgE), [ck(-2.6e-3,16,+1), ck(0,16,-1), ck(2.6e-3,16,+1)], geo);
A{end+1} = arrangement('a2', '(a′) 轴向环→← 相斥 d=10.4 缝内 13 匝', ...
    [mk('A',-d1/2,+1,mgE), mk('A',+d1/2,-1,mgE)], ck(0,13,+1), geo);
A{end+1} = arrangement('a3', '(a″) 轴向环→→ d=8.7 磁下（现状若为轴向）', ...
    [mk('A',-d0/2,+1,mgE), mk('A',+d0/2,+1,mgE)], [ck(-d0/2,29,+1), ck(+d0/2,29,+1)], geo);

%% 4. 双极化探头：D0（现有两径向环反极性）、D1（三径向环）、D2（三轴向环）、D1T（三厚径向环）
m2R = [mk('R',-d1/2,+1,mgE), mk('R',+d1/2,-1,mgE)];
m3R = [mk('R',-d1,+1,mgE), mk('R',0,-1,mgE), mk('R',+d1,+1,mgE)];
m3A = [mk('A',-d1,+1,mgE), mk('A',0,-1,mgE), mk('A',+d1,+1,mgE)];
mgT = struct('ri', 7.5e-3, 'ro', 12.5e-3, 'W', 8.0e-3);           % 拟购厚环：内径 15、外径 25、宽 8
m3T = [mk('R',-d1,+1,mgT), mk('R',0,-1,mgT), mk('R',+d1,+1,mgT)];
D = {};
D{end+1} = arrangement('D0w', 'D0 壁组 径向环+− 磁下 2×29 匝 (+−)', m2R, [ck(-d1/2,29,+1), ck(d1/2,29,-1)], geo);
D{end+1} = arrangement('D0l', 'D0 液组 径向环+− 缝内 13 匝', m2R, ck(0,13,+1), geo);
D{end+1} = arrangement('D1w', 'D1 壁组 径向环+−+ 磁下 3×29 匝 (+−+)', m3R, [ck(-d1,29,+1), ck(0,29,-1), ck(d1,29,+1)], geo);
D{end+1} = arrangement('D1l', 'D1 液组 径向环+−+ 缝内 2×13 匝 (+−)', m3R, [ck(-d1/2,13,+1), ck(d1/2,13,-1)], geo);
D{end+1} = arrangement('D2w', 'D2 壁组 轴向环→←→ 缝内 2×13 匝 (+−)', m3A, [ck(-d1/2,13,+1), ck(d1/2,13,-1)], geo);
D{end+1} = arrangement('D2l', 'D2 液组 轴向环→←→ 磁下 3×29 匝 (+−+)', m3A, [ck(-d1,29,+1), ck(0,29,-1), ck(d1,29,+1)], geo);
D{end+1} = arrangement('D1Tw', 'D1T 壁组 厚径向环(7.5/12.5)+−+ 磁下 3×29 匝', m3T, [ck(-d1,29,+1), ck(0,29,-1), ck(d1,29,+1)], geo);
D{end+1} = arrangement('D1Tl', 'D1T 液组 厚径向环(7.5/12.5)+−+ 缝内 2×13 匝', m3T, [ck(-d1/2,13,+1), ck(d1/2,13,-1)], geo);

%% 4b. 相邻磁体间的轴向力（给 E4 支架设计；正 = 相斥）
Fm.RR_same_87  = ringForceZ('R', mgE, -d0/2, +1, 'R', mgE, +d0/2, +1, geo.M);
Fm.RR_same_104 = ringForceZ('R', mgE, -d1/2, +1, 'R', mgE, +d1/2, +1, geo.M);
Fm.RR_opp_104  = ringForceZ('R', mgE, -d1/2, +1, 'R', mgE, +d1/2, -1, geo.M);
Fm.AA_rep_104  = ringForceZ('A', mgE, -d1/2, +1, 'A', mgE, +d1/2, -1, geo.M);
Fm.AA_att_104  = ringForceZ('A', mgE, -d1/2, +1, 'A', mgE, +d1/2, +1, geo.M);
Fm.TT_opp_104  = ringForceZ('R', mgT, -d1/2, +1, 'R', mgT, +d1/2, -1, geo.M);

%% 5. 充磁方向判别预测（高斯计）
rg = 7.5e-3;  rf = (mgE.ri + mgE.ro)/2;  zf = mgE.W/2 + 0.5e-3;
[BrR0, BzR0] = ringB('R', mgE, 0, +1, geo.M, rg, 0);
[BrRf, BzRf] = ringB('R', mgE, 0, +1, geo.M, rf, zf);
[BrA0, BzA0] = ringB('A', mgE, 0, +1, geo.M, rg, 0);
[BrAf, BzAf] = ringB('A', mgE, 0, +1, geo.M, rf, zf);
[~, BzAax]   = ringB('A', mgE, 0, +1, geo.M, 1e-9, 0);
zq = (-16:0.05:16)'*1e-3;
[BrRq, BzRq] = ringB('R', mgE, 0, +1, geo.M, rg, zq);
[BrAq, BzAq] = ringB('A', mgE, 0, +1, geo.M, rg, zq);

%% 6. 磁体几何变体
ri_list = [7.25 7.5 8.0]*1e-3;
t_list  = (2:1:7)*1e-3;
V.BrR = zeros(numel(ri_list), numel(t_list));  V.BzA = V.BrR;  V.BzRpair = V.BrR;  V.BrApair = V.BrR;
for i = 1:numel(ri_list)
    for j = 1:numel(t_list)
        mg = struct('ri', ri_list(i), 'ro', ri_list(i) + t_list(j), 'W', 8e-3);
        [br, ~]  = ringB('R', mg, 0, +1, geo.M, geo.rsurf, 0);            V.BrR(i,j) = abs(br);
        [~, bz]  = ringB('A', mg, 0, +1, geo.M, geo.rsurf, 0);            V.BzA(i,j) = abs(bz);
        [~, b1]  = ringB('R', mg, -d1/2, +1, geo.M, geo.rsurf, 0);
        [~, b2]  = ringB('R', mg, +d1/2, -1, geo.M, geo.rsurf, 0);        V.BzRpair(i,j) = abs(b1 + b2);
        [b1, ~]  = ringB('A', mg, -d1/2, +1, geo.M, geo.rsurf, 0);
        [b2, ~]  = ringB('A', mg, +d1/2, -1, geo.M, geo.rsurf, 0);        V.BrApair(i,j) = abs(b1 + b2);
    end
end
W_list = [6 8 10 12]*1e-3;  V.W_BrR = zeros(size(W_list));  V.W_BzA = V.W_BrR;
for j = 1:numel(W_list)
    mg = struct('ri', 8e-3, 'ro', 10e-3, 'W', W_list(j));
    [br, ~] = ringB('R', mg, 0, +1, geo.M, geo.rsurf, 0);  V.W_BrR(j) = abs(br);
    [~, bz] = ringB('A', mg, 0, +1, geo.M, geo.rsurf, 0);  V.W_BzA(j) = abs(bz);
end

%% 7. 图 E1-1 场分布
fig1 = figure('Visible', 'on', 'Name', 'E1-1 fields', 'Position', [40 40 1400 640], 'Color', 'w');
tl = tiledlayout(fig1, 2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:numel(A)
    ax = nexttile(tl);  plotFields(ax, A{i}, geo);
end
ax = nexttile(tl);
h1 = plot(ax, zq*1e3, BrRq, '-', 'Color', [0 0.3 0.8], 'LineWidth', 1.2); hold(ax, 'on');
h2 = plot(ax, zq*1e3, BzRq, '--', 'Color', [0 0.3 0.8], 'LineWidth', 1.2);
h3 = plot(ax, zq*1e3, BrAq, '-', 'Color', [0.85 0.2 0.1], 'LineWidth', 1.2);
h4 = plot(ax, zq*1e3, BzAq, '--', 'Color', [0.85 0.2 0.1], 'LineWidth', 1.2);
xlim(ax, [-16 16]); xlabel(ax, 'z / mm'); ylabel(ax, 'B / T');
legend(ax, [h1 h2 h3 h4], {'径向充磁 B_r', '径向充磁 B_z', '轴向充磁 B_r', '轴向充磁 B_z'}, ...
    'Location', 'northeast', 'Box', 'off', 'FontSize', 8);
title(ax, '高斯计判别：内孔 r = 7.5 mm 处（单环）', 'FontWeight', 'normal', 'FontSize', 9);
set(ax, 'Box', 'on', 'Layer', 'top');
hideToolbars(fig1);
exportgraphics(fig1, fullfile(figdir, 'E1_1_fields.png'), 'Resolution', 200);
savefig(fig1, fullfile(figdir, 'E1_1_fields.fig'));

%% 8. 图 E1-2 力分布与空间谱
sel = {'a1', 'b1', 'c2', 'a2'};
fig2 = figure('Visible', 'on', 'Name', 'E1-2 forces', 'Position', [60 60 1400 640], 'Color', 'w');
tl = tiledlayout(fig2, 2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:numel(sel)
    S = A{strcmp(cellfun(@(s) s.id, A, 'UniformOutput', false), sel{i})};
    ax = nexttile(tl, i);      plotForce(ax, S, geo);
    ax = nexttile(tl, i + 4);  plotSpec(ax, S, geo);
end
hideToolbars(fig2);
exportgraphics(fig2, fullfile(figdir, 'E1_2_forces.png'), 'Resolution', 200);
savefig(fig2, fullfile(figdir, 'E1_2_forces.fig'));

%% 9. 图 E1-3 双极化探头 D0 / D1 / D2
fig3 = figure('Visible', 'on', 'Name', 'E1-3 dual', 'Position', [80 80 1400 900], 'Color', 'w');
tl = tiledlayout(fig3, 3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
pairs = {D{1}, D{2}; D{3}, D{4}; D{5}, D{6}};
for c = 1:3
    Sw = pairs{c,1}; Sl = pairs{c,2};
    ax = nexttile(tl, c);  plotFields(ax, Sw, geo);
    ax = nexttile(tl, c + 3);
    z = geo.z*1e3;
    h1 = plot(ax, z, Sw.fz, '-',  'Color', [0 0.3 0.8],    'LineWidth', 1.2); hold(ax, 'on');
    h2 = plot(ax, z, Sw.fr, '--', 'Color', [0 0.3 0.8],    'LineWidth', 1.0);
    h3 = plot(ax, z, Sl.fr, '-',  'Color', [0.85 0.2 0.1], 'LineWidth', 1.2);
    h4 = plot(ax, z, Sl.fz, '--', 'Color', [0.85 0.2 0.1], 'LineWidth', 1.0);
    yl = max(abs([Sw.fz; Sw.fr; Sl.fr; Sl.fz]))*1.25;  ylim(ax, [-yl yl]); xlim(ax, [-25 25]);
    xlabel(ax, 'z / mm'); ylabel(ax, 'f / (N·m^{-1}·A^{-1})');
    legend(ax, [h1 h2 h3 h4], {'壁组 f_z', '壁组 f_r（串扰）', '液组 f_r', '液组 f_z（串扰）'}, ...
        'Location', 'southwest', 'Box', 'off', 'FontSize', 8);
    title(ax, [Sw.id(1:2) ' 两组线圈的力分布'], 'FontWeight', 'normal', 'FontSize', 9);
    set(ax, 'Box', 'on', 'Layer', 'top');
    ax = nexttile(tl, c + 6);
    lam = geo.lam*1e3;
    h1 = plot(ax, lam, Sw.Fz, '-',  'Color', [0 0.3 0.8],    'LineWidth', 1.2); hold(ax, 'on');
    h2 = plot(ax, lam, Sw.Fr, '--', 'Color', [0 0.3 0.8],    'LineWidth', 1.0);
    h3 = plot(ax, lam, Sl.Fr, '-',  'Color', [0.85 0.2 0.1], 'LineWidth', 1.2);
    h4 = plot(ax, lam, Sl.Fz, '--', 'Color', [0.85 0.2 0.1], 'LineWidth', 1.0);
    yl = max([Sw.Fz Sw.Fr Sl.Fr Sl.Fz])*1.15;  ylim(ax, [0 yl]); xlim(ax, [2 30]);
    markTargets(ax, geo);
    xlabel(ax, 'λ / mm'); ylabel(ax, '|F(λ)| / (N·A^{-1})');
    legend(ax, [h1 h2 h3 h4], {'壁组 F_z', '壁组 F_r', '液组 F_r', '液组 F_z'}, 'Location', 'northeast', 'Box', 'off', 'FontSize', 8);
    title(ax, [Sw.id(1:2) ' 空间谱'], 'FontWeight', 'normal', 'FontSize', 9);
    set(ax, 'Box', 'on', 'Layer', 'top');
end
hideToolbars(fig3);
exportgraphics(fig3, fullfile(figdir, 'E1_3_dual.png'), 'Resolution', 200);
savefig(fig3, fullfile(figdir, 'E1_3_dual.fig'));

%% 10. 图 E1-4 磁体几何变体
fig4 = figure('Visible', 'on', 'Name', 'E1-4 variants', 'Position', [100 100 1400 380], 'Color', 'w');
tl = tiledlayout(fig4, 1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
names = {'径向环 磁下 |B_r|', '轴向环 磁下 |B_z|', '径向环+− 缝内 |B_z|', '轴向环→← 缝内 |B_r|'};
dat   = {V.BrR, V.BzA, V.BzRpair, V.BrApair};
cols  = [0 0.3 0.8; 0.85 0.2 0.1; 0.1 0.6 0.3];
for k = 1:4
    ax = nexttile(tl, k);  hold(ax, 'on');
    for i = 1:numel(ri_list)
        plot(ax, t_list*1e3, dat{k}(i,:), '-o', 'Color', cols(i,:), 'LineWidth', 1.2, 'MarkerSize', 4, 'MarkerFaceColor', cols(i,:));
    end
    plot(ax, 2, dat{k}(3,1), 'ks', 'MarkerSize', 9, 'LineWidth', 1.2);
    xlabel(ax, '径向厚度 t_m / mm'); ylabel(ax, 'B 在管外壁 / T');
    title(ax, names{k}, 'FontWeight', 'normal', 'FontSize', 9);
    legend(ax, {'内半径 7.25', '内半径 7.5', '内半径 8.0（现有）', '现有磁体'}, 'Location', 'southeast', 'Box', 'off', 'FontSize', 8);
    set(ax, 'Box', 'on', 'Layer', 'top'); xlim(ax, [1.5 7.5]);
end
hideToolbars(fig4);
exportgraphics(fig4, fullfile(figdir, 'E1_4_variants.png'), 'Resolution', 200);
savefig(fig4, fullfile(figdir, 'E1_4_variants.fig'));

%% 11. 结果文本
fid = fopen(fullfile(here, 'E1_结果.txt'), 'w', 'n', 'UTF-8');
fprintf(fid, 'E1 力场方案 结果  %s\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')));
fprintf(fid, '输入：磁体 N45H B_rem = %.2f T，现有环 内半径 %.2f / 外半径 %.2f / 宽 %.1f mm；管外壁 r_o = %.1f mm；\n', ...
    Brem, mgE.ri*1e3, mgE.ro*1e3, mgE.W*1e3, geo.rsurf*1e3);
fprintf(fid, '      线圈 0.15 mm 单层，导线中心离壁 %.1f mm；涡流镜像近似；力 = K × B，每安培线圈电流；周长 2π r_o 已乘。\n', geo.h*1e3);
fprintf(fid, '自检：单环轴上场相对误差 %.1e；轴向环轴上场相对误差 %.1e\n\n', err1, err2);
fprintf(fid, '== 现有磁体（8×2 mm，内半径 8）七种布置：线圈区内 |B| 峰值 / 空间谱 F(λ) N/A / 纯度 dB ==\n');
fprintf(fid, '%-6s %-44s %7s %7s | %8s %8s %6s | %8s %8s %6s | %8s %8s %6s | %8s %8s %6s | %7s %7s %5s\n', ...
    'id', '布置', 'Br_coil', 'Bz_coil', 'Fz10.4', 'Fr10.4', 'ax dB', 'Fz6.3', 'Fr6.3', 'rd dB', ...
    'Fz5.2', 'Fr5.2', 'rd dB', 'Fz4.5', 'Fr4.5', 'rd dB', 'int|fz|', 'int|fr|', 'N');
for i = 1:numel(A), writeRow(fid, A{i}); end
fprintf(fid, '\n== 双极化探头（节距 10.4 mm；D0 两环、D1/D2 三环、D1T 三厚环 内径 15/外径 25/宽 8）==\n');
for i = 1:numel(D), writeRow(fid, D{i}); end
fprintf(fid, '\n纯度定义：ax dB = 20 log10(Fz/Fr)（轴向力用，越大越纯）；rd dB = 20 log10(Fr/Fz)（径向力用）。\n');
fprintf(fid, 'Fλ 为 |∫ f e^{−ikz} dz|，单位 N 每安培线圈电流；int|f| 为 ∫|f| dz（N/A）；N 为总匝数。\n');
fprintf(fid, '\n== 相邻磁体间轴向力（N，正 = 相斥，负 = 相吸；给支架设计）==\n');
fprintf(fid, '径向环 同极性 d = 8.7（现状）：%+.1f；同极性 d = 10.4：%+.1f；反极性 d = 10.4：%+.1f\n', Fm.RR_same_87, Fm.RR_same_104, Fm.RR_opp_104);
fprintf(fid, '轴向环 相斥(→←) d = 10.4：%+.1f；相吸(→→) d = 10.4：%+.1f；厚径向环 反极性 d = 10.4：%+.1f\n', Fm.AA_rep_104, Fm.AA_att_104, Fm.TT_opp_104);
fprintf(fid, '\n== 充磁方向判别预测（B_rem = 1.35 T；高斯计探头离表面 0.5 mm）==\n');
fprintf(fid, '径向充磁：内孔中面 r = 7.5 mm：B_r = %+.3f T（法向），B_z = %+.3f T；端面外 0.5 mm（r = 9）：B_z = %+.3f T（法向），B_r = %+.3f T\n', BrR0, BzR0, BzRf, BrRf);
fprintf(fid, '轴向充磁：内孔中面 r = 7.5 mm：B_r = %+.3f T，B_z = %+.3f T（切向，与磁化反向）；端面外 0.5 mm：B_z = %+.3f T（法向），B_r = %+.3f T；孔中心轴上 B_z = %+.3f T\n', BrA0, BzA0, BzAf, BrAf, BzAax);
fprintf(fid, '直径充磁（若买到的是这种）：内孔场随周向角按 cos φ 变化，转一圈变号一次；径向充磁转一圈读数不变。\n');
fprintf(fid, '\n== 磁体几何变体：管外壁处 |B| (T)，W = 8 mm，列 = 径向厚度 %s mm ==\n', mat2str(t_list*1e3));
lab = {'径向环 磁下 |B_r|', '轴向环 磁下 |B_z|', '径向环+− d=10.4 缝内 |B_z|', '轴向环→← d=10.4 缝内 |B_r|'};
for k = 1:4
    fprintf(fid, '%s\n', lab{k});
    for i = 1:numel(ri_list)
        fprintf(fid, '  内半径 %.2f：%s\n', ri_list(i)*1e3, sprintf(' %.3f', dat{k}(i,:)));
    end
end
fprintf(fid, '宽度 W 变体（内半径 8、厚 2）：W = %s mm → 径向环 |B_r| %s；轴向环 |B_z| %s\n', ...
    mat2str(W_list*1e3), sprintf(' %.3f', V.W_BrR), sprintf(' %.3f', V.W_BzA));
fclose(fid);

%% 12. 给 E2/E3 的场数据
fields = struct('z', geo.z, 'rsurf', geo.rsurf, 'Brem', Brem, 'mgE', mgE, 'geo', geo);
for i = 1:numel(A), fields.(A{i}.id) = rmfield(A{i}, {'Fz', 'Fr'}); end
for i = 1:numel(D), fields.(D{i}.id) = rmfield(D{i}, {'Fz', 'Fr'}); end
fields.variants = V;  fields.ri_list = ri_list;  fields.t_list = t_list;  fields.W_list = W_list;
fields.magnetForce = Fm;  fields.mgT = mgT;
save(fullfile(here, 'E1_fields.mat'), '-struct', 'fields');

figure(fig1); figure(fig2); figure(fig3); figure(fig4); drawnow;
fprintf('E1 完成：图与结果已写入 %s\n', here);
end

%% ---------- 局部函数 ----------
function S = arrangement(id, name, mags, coils, geo)
S.id = id;  S.name = name;  S.mags = mags;  S.coils = coils;
Br = 0;  Bz = 0;
for i = 1:numel(mags)
    [br, bz] = ringB(mags(i).type, mags(i).mg, mags(i).zc, mags(i).sgn, geo.M, geo.rsurf, geo.z);
    Br = Br + br;  Bz = Bz + bz;
end
K = 0;
for i = 1:numel(coils)
    K = K + coilK(coils(i).zc, coils(i).Nt, coils(i).sgn, geo.rc, geo.rsurf, geo.pitch, geo.z);
end
S.Br = Br;  S.Bz = Bz;  S.K = K;
S.fz = -2*pi*geo.rsurf*K.*Br;                 % N/m 每 A
S.fr =  2*pi*geo.rsurf*K.*Bz;
S.Fz   = spec(S.fz, geo.z, geo.lam);    S.Fr   = spec(S.fr, geo.z, geo.lam);
S.Fz_t = spec(S.fz, geo.z, geo.lam_t);  S.Fr_t = spec(S.fr, geo.z, geo.lam_t);
S.Iz = trapz(geo.z, abs(S.fz));  S.Ir = trapz(geo.z, abs(S.fr));
mask = false(size(geo.z));
for i = 1:numel(coils)
    mask = mask | (abs(geo.z - coils(i).zc) <= coils(i).Nt*geo.pitch/2);
end
S.Br_coil = max(abs(Br(mask)));  S.Bz_coil = max(abs(Bz(mask)));
S.Ntot = sum([coils.Nt]);
end

function F = spec(f, z, lam)
k = 2*pi./lam(:).';
F = abs(trapz(z, f.*exp(-1i*z*k), 1));
end

function Fz = ringForceZ(type1, mg1, zc1, sgn1, type2, mg2, zc2, sgn2, M)
% 环 2 受环 1 的轴向力：环 2 的等效电流环在环 1 场中受力 F_z = −2π a I B_r,1 叠加（N）。
% 返回值按"环 1 在左（zc1 < zc2）"解释：正 = 相斥，负 = 相吸。
Fz = 0;
switch upper(type2)
    case 'R'
        na = max(40, round((mg2.ro - mg2.ri)/0.05e-3));
        ae = linspace(mg2.ri, mg2.ro, na+1);  a = (ae(1:end-1) + ae(2:end))/2;  da = (mg2.ro - mg2.ri)/na;
        for k = 1:na
            [br, ~] = ringB(type1, mg1, zc1, sgn1, M, a(k), zc2 - mg2.W/2);  Fz = Fz - 2*pi*a(k)*(+sgn2*M*da)*br;
            [br, ~] = ringB(type1, mg1, zc1, sgn1, M, a(k), zc2 + mg2.W/2);  Fz = Fz - 2*pi*a(k)*(-sgn2*M*da)*br;
        end
    case 'A'
        nz = max(80, round(mg2.W/0.1e-3));
        ze = linspace(zc2 - mg2.W/2, zc2 + mg2.W/2, nz+1);  zs = (ze(1:end-1) + ze(2:end))/2;  dz = mg2.W/nz;
        for k = 1:nz
            [br, ~] = ringB(type1, mg1, zc1, sgn1, M, mg2.ro, zs(k));  Fz = Fz - 2*pi*mg2.ro*(+sgn2*M*dz)*br;
            [br, ~] = ringB(type1, mg1, zc1, sgn1, M, mg2.ri, zs(k));  Fz = Fz - 2*pi*mg2.ri*(-sgn2*M*dz)*br;
        end
end
end

function hideToolbars(fig)
axs = findall(fig, 'Type', 'axes');
for k = 1:numel(axs)
    try, axs(k).Toolbar.Visible = 'off'; catch, end
end
end

function writeRow(fid, S)
pur = @(num, den) 20*log10(max(num, 1e-12)/max(den, 1e-12));
t = S.Fz_t;  r = S.Fr_t;
fprintf(fid, '%-6s %-44s %7.3f %7.3f | %8.4f %8.4f %6.1f | %8.4f %8.4f %6.1f | %8.4f %8.4f %6.1f | %8.4f %8.4f %6.1f | %7.3f %7.3f %5d\n', ...
    S.id, S.name, S.Br_coil, S.Bz_coil, t(1), r(1), pur(t(1), r(1)), t(2), r(2), pur(r(2), t(2)), ...
    t(3), r(3), pur(r(3), t(3)), t(4), r(4), pur(r(4), t(4)), S.Iz, S.Ir, S.Ntot);
end

function drawLayout(ax, S, geo, yl)
for i = 1:numel(S.mags)
    m = S.mags(i);  x0 = (m.zc - m.mg.W/2)*1e3;  w = m.mg.W*1e3;
    p = patch(ax, [x0 x0+w x0+w x0], [yl(1) yl(1) yl(2) yl(2)], [0.88 0.88 0.88], ...
        'EdgeColor', [0.6 0.6 0.6], 'FaceAlpha', 0.6);
    uistack(p, 'bottom');
    if upper(m.type) == 'R'
        if m.sgn > 0, lab = 'R+'; else, lab = 'R−'; end
    else
        if m.sgn > 0, lab = 'A→'; else, lab = 'A←'; end
    end
    text(ax, m.zc*1e3, yl(2)*0.88, lab, 'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', [0.3 0.3 0.3]);
end
for i = 1:numel(S.coils)
    c = S.coils(i);  w = c.Nt*geo.pitch*1e3;
    plot(ax, [c.zc*1e3 - w/2, c.zc*1e3 + w/2], [yl(1) yl(1)]*0.9, 'k-', 'LineWidth', 3);
    if c.sgn > 0, lab = '+'; else, lab = '−'; end
    text(ax, c.zc*1e3, yl(1)*0.76, lab, 'HorizontalAlignment', 'center', 'FontSize', 9);
end
end

function plotFields(ax, S, geo)
z = geo.z*1e3;
h1 = plot(ax, z, S.Br, '-', 'Color', [0 0.3 0.8], 'LineWidth', 1.2);  hold(ax, 'on');
h2 = plot(ax, z, S.Bz, '-', 'Color', [0.85 0.2 0.1], 'LineWidth', 1.2);
yl = max([abs(S.Br); abs(S.Bz)])*1.25;  if yl == 0, yl = 0.1; end
ylim(ax, [-yl yl]);  xlim(ax, [-25 25]);
drawLayout(ax, S, geo, [-yl yl]);
xlabel(ax, 'z / mm');  ylabel(ax, 'B 在管外壁 / T');
legend(ax, [h1 h2], {'B_r', 'B_z'}, 'Location', 'northeast', 'Box', 'off', 'FontSize', 8);
title(ax, S.name, 'FontWeight', 'normal', 'FontSize', 9);
set(ax, 'Box', 'on', 'Layer', 'top');
end

function plotForce(ax, S, geo)
z = geo.z*1e3;
h1 = plot(ax, z, S.fz, '-', 'Color', [0 0.3 0.8], 'LineWidth', 1.2);  hold(ax, 'on');
h2 = plot(ax, z, S.fr, '-', 'Color', [0.85 0.2 0.1], 'LineWidth', 1.2);
yl = max([abs(S.fz); abs(S.fr)])*1.25;  if yl == 0, yl = 1; end
ylim(ax, [-yl yl]);  xlim(ax, [-25 25]);
drawLayout(ax, S, geo, [-yl yl]);
xlabel(ax, 'z / mm');  ylabel(ax, 'f / (N·m^{-1}·A^{-1})');
legend(ax, [h1 h2], {'f_z 轴向', 'f_r 径向'}, 'Location', 'northeast', 'Box', 'off', 'FontSize', 8);
title(ax, S.name, 'FontWeight', 'normal', 'FontSize', 9);
set(ax, 'Box', 'on', 'Layer', 'top');
end

function plotSpec(ax, S, geo)
lam = geo.lam*1e3;
plot(ax, lam, S.Fz, '-', 'Color', [0 0.3 0.8], 'LineWidth', 1.2);  hold(ax, 'on');
plot(ax, lam, S.Fr, '-', 'Color', [0.85 0.2 0.1], 'LineWidth', 1.2);
yl = max([S.Fz S.Fr])*1.15;  if yl == 0, yl = 1; end
ylim(ax, [0 yl]);  xlim(ax, [2 30]);
markTargets(ax, geo);
xlabel(ax, 'λ / mm');  ylabel(ax, '|F(λ)| / (N·A^{-1})');
legend(ax, {'F_z 轴向', 'F_r 径向'}, 'Location', 'northeast', 'Box', 'off', 'FontSize', 8);
set(ax, 'Box', 'on', 'Layer', 'top');
end

function markTargets(ax, geo)
for j = 1:numel(geo.lam_t)
    xline(ax, geo.lam_t(j)*1e3, ':', sprintf('%.1f', geo.lam_t(j)*1e3), 'Color', [0.4 0.4 0.4], ...
        'LabelOrientation', 'aligned', 'LabelVerticalAlignment', 'top', ...
        'LabelHorizontalAlignment', 'left', 'FontSize', 8);
end
end
