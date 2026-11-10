function E2_coil_spectrum()
%E2_COIL_SPECTRUM  E2 线圈几何：元件宽、节距、匝数的空间谱优化（用 E1 的磁场）
% 输入来源：
%   E1_fields.mat（E1_force_field.m 输出）：D0 / D1 / D1T 三种磁体叠列在管外壁 r_o = 6.5 mm 处的 B_r(z)、B_z(z)，
%       评估轴 z、线圈几何 geo（导线中心半径 7.0 mm、匝距 0.15 mm）、目标波长 [10.4 6.3 5.2 4.5] mm
%   coilK.m：周向线圈在管壁感应的面电流（镜像近似，每安培）
% 方法：固定磁体叠列，扫描壁组元件宽 w_w（磁下，电流随磁体极性交替）与液组元件宽 w_l（缝下，两缝电流反向），
%   匝数 N = round(w/0.15) 单层；力 f_z = −2π r_o K B_r、f_r = 2π r_o K B_z；空间谱 F(λ) = |∫ f e^{−i2πz/λ} dz|。
%   按互易性，F(λ) 同时是发射力幅（N/A）和接收灵敏度（V 每单位模态速度幅）的度量，所以直接用 F 选宽度。
% 符号与单位：同 E1；mm、N/A、µH。
% 输出：figures\E2_1_coil_sweep.*、E2_2_final_spectra.*；E2_结果.txt；E2_coils.mat（最终线圈规格，给 E3/E4）
% 运行：在本目录 matlab -batch "E2_coil_spectrum"，约 1 min。

here   = fileparts(mfilename('fullpath'));
figdir = fullfile(here, 'figures');
addpath(here);
F   = load(fullfile(here, 'E1_fields.mat'));
geo = F.geo;  z = F.z;  lam = geo.lam;  lam_t = geo.lam_t;  pitch = geo.pitch;  r_o = geo.rsurf;

FONT = 'TimesSimSun';
try
    if any(strcmpi(FONT, listfonts))
        set(groot, 'defaultAxesFontName', FONT, 'defaultTextFontName', FONT, 'defaultLegendFontName', FONT);
    end
catch
end
set(groot, 'defaultTextInterpreter', 'tex', 'defaultAxesTickLabelInterpreter', 'tex', ...
    'defaultLegendInterpreter', 'tex', 'defaultAxesFontSize', 9);

%% 1. 三种磁体叠列（场来自 E1）与线圈位置
d1 = 10.4e-3;
St{1} = struct('id', 'D0',  'Br', F.D0w.Br,  'Bz', F.D0w.Bz,  'zw', [-d1/2 d1/2], 'sw', [1 -1],   'zl', 0,            'sl', 1);
St{2} = struct('id', 'D1',  'Br', F.D1w.Br,  'Bz', F.D1w.Bz,  'zw', [-d1 0 d1],   'sw', [1 -1 1], 'zl', [-d1/2 d1/2], 'sl', [1 -1]);
St{3} = struct('id', 'D1T', 'Br', F.D1Tw.Br, 'Bz', F.D1Tw.Bz, 'zw', [-d1 0 d1],   'sw', [1 -1 1], 'zl', [-d1/2 d1/2], 'sl', [1 -1]);

%% 2. 扫描元件宽
ww = (1.5:0.15:7.05)*1e-3;      % 壁组元件宽
wl = (0.9:0.15:4.05)*1e-3;      % 液组元件宽
for s = 1:numel(St)
    S = St{s};
    S.Fw = zeros(numel(ww), 4);  S.Fw_x = S.Fw;  S.Nw = zeros(size(ww));
    for i = 1:numel(ww)
        N = max(1, round(ww(i)/pitch));  S.Nw(i) = N;
        K = 0;  for e = 1:numel(S.zw), K = K + coilK(S.zw(e), N, S.sw(e), geo.rc, r_o, pitch, z); end
        S.Fw(i,:)   = spec(-2*pi*r_o*K.*S.Br, z, lam_t);      % 轴向（所需）
        S.Fw_x(i,:) = spec( 2*pi*r_o*K.*S.Bz, z, lam_t);      % 径向（反串）
    end
    S.Fl = zeros(numel(wl), 4);  S.Fl_x = S.Fl;  S.Nl = zeros(size(wl));
    for i = 1:numel(wl)
        N = max(1, round(wl(i)/pitch));  S.Nl(i) = N;
        K = 0;  for e = 1:numel(S.zl), K = K + coilK(S.zl(e), N, S.sl(e), geo.rc, r_o, pitch, z); end
        S.Fl(i,:)   = spec( 2*pi*r_o*K.*S.Bz, z, lam_t);      % 径向（所需）
        S.Fl_x(i,:) = spec(-2*pi*r_o*K.*S.Br, z, lam_t);      % 轴向（反串）
    end
    St{s} = S;
end

%% 3. 选定宽度：壁组取 F_z(10.4) 最大；液组取 F_r(5.2) 最大（同时看 6.3 / 4.5 不能塌）
for s = 1:numel(St)
    [~, iw] = max(St{s}.Fw(:,1));      St{s}.ww_opt = ww(iw);  St{s}.Nw_opt = St{s}.Nw(iw);
    [~, il] = max(St{s}.Fl(:,3));      St{s}.wl_opt = wl(il);  St{s}.Nl_opt = St{s}.Nl(il);
    St{s}.iw = iw;  St{s}.il = il;
end
% 统一规格：三种叠列取同一宽度（便于一套骨架）。
%   壁组不取扫描最大值（6.15 mm）：从 4.35 到 6.15 mm F_z(10.4) 只涨 9%，径向反串翻倍、匝数多 40%；
%   取经典 λ_w/2 = 5.2 mm（35 匝），F_z 达最大值的 98%，且不出 B_r 平台（±3 mm）。
%   液组取 F_r(5.2) 最大值（2.85 mm，19 匝），该极值很平（2.4–3.3 mm 内差 <2%）。
ww_f = 5.2e-3;         Nw_f = round(ww_f/pitch);   ww_f = Nw_f*pitch;
wl_f = St{3}.wl_opt;   Nl_f = round(wl_f/pitch);
[~, iwf] = min(abs(ww - ww_f));
ratio_w = [St{1}.Fw(iwf,1)/max(St{1}.Fw(:,1)), St{2}.Fw(iwf,1)/max(St{2}.Fw(:,1)), St{3}.Fw(iwf,1)/max(St{3}.Fw(:,1))];
% 壁模态径向反串的作用（薄板近似）：S0 外表面 u_r = −[ν h/(2(1−ν))] ∂u_z/∂z，与 u_z 正交相位；
%   叠列对称（f_z 偶、f_r 奇）→ 两个传播方向幅值相同；反串项与轴向项同号相加，增益 ≈ 1 + η·F_r/F_z
nu = 0.33;  hwall = 1e-3;  kw = 2*pi/10.4e-3;
eta_plate = nu*hwall*kw/(2*(1 - nu));

%% 4. 最终设计的完整谱（两组、两分量）+ "一套线圈兼顾两族"的检验
for s = 1:numel(St)
    S = St{s};
    Kw = 0;  for e = 1:numel(S.zw), Kw = Kw + coilK(S.zw(e), Nw_f, S.sw(e), geo.rc, r_o, pitch, z); end
    Kl = 0;  for e = 1:numel(S.zl), Kl = Kl + coilK(S.zl(e), Nl_f, S.sl(e), geo.rc, r_o, pitch, z); end
    S.fzw = -2*pi*r_o*Kw.*S.Br;  S.frw = 2*pi*r_o*Kw.*S.Bz;
    S.frl =  2*pi*r_o*Kl.*S.Bz;  S.fzl = -2*pi*r_o*Kl.*S.Br;
    S.Fzw = spec(S.fzw, z, lam);  S.Frw = spec(S.frw, z, lam);
    S.Frl = spec(S.frl, z, lam);  S.Fzl = spec(S.fzl, z, lam);
    S.Fzw_t = spec(S.fzw, z, lam_t);  S.Frw_t = spec(S.frw, z, lam_t);
    S.Frl_t = spec(S.frl, z, lam_t);  S.Fzl_t = spec(S.fzl, z, lam_t);
    % 单套方案：磁下窄元件（2.6 mm，17 匝），看它在 λ = 5.2 的轴向力有多大（方向仍是轴向）
    Ks = 0;  for e = 1:numel(S.zw), Ks = Ks + coilK(S.zw(e), 17, S.sw(e), geo.rc, r_o, pitch, z); end
    S.Fz_single_t = spec(-2*pi*r_o*Ks.*S.Br, z, lam_t);
    S.Fr_single_t = spec( 2*pi*r_o*Ks.*S.Bz, z, lam_t);
    St{s} = S;
end

%% 5. 电感粗估（Wheeler 单层空心螺线管；管壁涡流镜像会再降 30–50%，以 LCR 表实测为准）
rin = geo.rc/0.0254;
Lw = @(N, w) (rin^2*N^2)/(9*rin + 10*w/0.0254);        % µH
L_w1 = Lw(Nw_f, ww_f);  L_l1 = Lw(Nl_f, wl_f);
Rdc = @(N) N*2*pi*geo.rc*1.72e-8/(pi*(0.075e-3)^2);      % 0.15 mm 铜线直流电阻 Ω

%% 6. 图 E2-1 扫描
fig1 = figure('Visible', 'on', 'Name', 'E2-1 coil sweep', 'Position', [40 40 1400 640], 'Color', 'w');
tl = tiledlayout(fig1, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
cols = [0.1 0.6 0.3; 0 0.3 0.8; 0.85 0.2 0.1];
ax = nexttile(tl, 1); hold(ax, 'on');
for s = 1:3
    plot(ax, ww*1e3, St{s}.Fw(:,1), '-', 'Color', cols(s,:), 'LineWidth', 1.3);
    plot(ax, ww*1e3, St{s}.Fw_x(:,1), '--', 'Color', cols(s,:), 'LineWidth', 1.0);
end
xline(ax, 4.35, ':', '现有 4.35', 'Color', [0.4 0.4 0.4], 'LabelOrientation', 'aligned', 'FontSize', 8);
xline(ax, ww_f*1e3, '-', sprintf('选 %.2f', ww_f*1e3), 'Color', [0.4 0.4 0.4], 'LabelOrientation', 'aligned', 'FontSize', 8);
xlabel(ax, '壁组元件宽 w_w / mm'); ylabel(ax, 'F(λ = 10.4) / (N·A^{-1})');
legend(ax, {'D0 F_z', 'D0 F_r 反串', 'D1 F_z', 'D1 F_r 反串', 'D1T F_z', 'D1T F_r 反串'}, 'Location', 'northwest', 'Box', 'off', 'FontSize', 8);
title(ax, '壁组：λ = 10.4 mm 的响应随元件宽', 'FontWeight', 'normal', 'FontSize', 9);
set(ax, 'Box', 'on', 'Layer', 'top');
ax = nexttile(tl, 2); hold(ax, 'on');
for s = 1:3
    plot(ax, ww*1e3, St{s}.Fw(:,1)./St{s}.Nw(:), '-', 'Color', cols(s,:), 'LineWidth', 1.3);
end
xlabel(ax, '壁组元件宽 w_w / mm'); ylabel(ax, 'F_z(10.4) 每匝 / (N·A^{-1})');
legend(ax, {'D0', 'D1', 'D1T'}, 'Location', 'northeast', 'Box', 'off', 'FontSize', 8);
title(ax, '壁组：每匝效率', 'FontWeight', 'normal', 'FontSize', 9);
set(ax, 'Box', 'on', 'Layer', 'top');
ax = nexttile(tl, 3); hold(ax, 'on');
s = 3;
plot(ax, ww*1e3, St{s}.Fw(:,1), '-', 'Color', cols(s,:), 'LineWidth', 1.3);
plot(ax, ww*1e3, St{s}.Fw(:,3), '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0);
plot(ax, ww*1e3, St{s}.Fw_x(:,3), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0);
xlabel(ax, '壁组元件宽 w_w / mm'); ylabel(ax, 'F / (N·A^{-1})');
legend(ax, {'F_z(10.4) 所需', 'F_z(5.2) 旁瓣', 'F_r(5.2) 旁瓣'}, 'Location', 'northwest', 'Box', 'off', 'FontSize', 8);
title(ax, 'D1T 壁组：λ = 5.2 处的旁瓣', 'FontWeight', 'normal', 'FontSize', 9);
set(ax, 'Box', 'on', 'Layer', 'top');
ax = nexttile(tl, 4); hold(ax, 'on');
for s = 1:3
    plot(ax, wl*1e3, St{s}.Fl(:,3), '-', 'Color', cols(s,:), 'LineWidth', 1.3);
    plot(ax, wl*1e3, St{s}.Fl_x(:,3), '--', 'Color', cols(s,:), 'LineWidth', 1.0);
end
xline(ax, wl_f*1e3, '-', sprintf('选 %.2f', wl_f*1e3), 'Color', [0.4 0.4 0.4], 'LabelOrientation', 'aligned', 'FontSize', 8);
xlabel(ax, '液组元件宽 w_l / mm'); ylabel(ax, 'F(λ = 5.2) / (N·A^{-1})');
legend(ax, {'D0 F_r', 'D0 F_z 反串', 'D1 F_r', 'D1 F_z 反串', 'D1T F_r', 'D1T F_z 反串'}, 'Location', 'northwest', 'Box', 'off', 'FontSize', 8);
title(ax, '液组：λ = 5.2 mm 的响应随元件宽', 'FontWeight', 'normal', 'FontSize', 9);
set(ax, 'Box', 'on', 'Layer', 'top');
ax = nexttile(tl, 5); hold(ax, 'on');
s = 3;
plot(ax, wl*1e3, St{s}.Fl(:,2), '-', 'Color', [0.1 0.6 0.3], 'LineWidth', 1.3);
plot(ax, wl*1e3, St{s}.Fl(:,3), '-', 'Color', [0.85 0.2 0.1], 'LineWidth', 1.3);
plot(ax, wl*1e3, St{s}.Fl(:,4), '-', 'Color', [0 0.3 0.8], 'LineWidth', 1.3);
xlabel(ax, '液组元件宽 w_l / mm'); ylabel(ax, 'F_r / (N·A^{-1})');
legend(ax, {'λ = 6.3（1904 m/s）', 'λ = 5.2（1560 m/s）', 'λ = 4.5（1359 m/s）'}, 'Location', 'northwest', 'Box', 'off', 'FontSize', 8);
title(ax, 'D1T 液组：300 kHz 三个液模态', 'FontWeight', 'normal', 'FontSize', 9);
set(ax, 'Box', 'on', 'Layer', 'top');
ax = nexttile(tl, 6); hold(ax, 'on');
for s = 1:3
    plot(ax, wl*1e3, St{s}.Fl(:,3)./St{s}.Nl(:), '-', 'Color', cols(s,:), 'LineWidth', 1.3);
end
xlabel(ax, '液组元件宽 w_l / mm'); ylabel(ax, 'F_r(5.2) 每匝 / (N·A^{-1})');
legend(ax, {'D0', 'D1', 'D1T'}, 'Location', 'northeast', 'Box', 'off', 'FontSize', 8);
title(ax, '液组：每匝效率', 'FontWeight', 'normal', 'FontSize', 9);
set(ax, 'Box', 'on', 'Layer', 'top');
hideToolbars(fig1);
exportgraphics(fig1, fullfile(figdir, 'E2_1_coil_sweep.png'), 'Resolution', 200);
savefig(fig1, fullfile(figdir, 'E2_1_coil_sweep.fig'));

%% 7. 图 E2-2 最终设计的空间谱（D0 与 D1T）
fig2 = figure('Visible', 'on', 'Name', 'E2-2 final spectra', 'Position', [60 60 1400 640], 'Color', 'w');
tl = tiledlayout(fig2, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for c = 1:2
    s = [1 3];  s = s(c);  S = St{s};
    ax = nexttile(tl, c); hold(ax, 'on');
    zz = z*1e3;
    h1 = plot(ax, zz, S.fzw, '-',  'Color', [0 0.3 0.8],    'LineWidth', 1.2);
    h2 = plot(ax, zz, S.frw, '--', 'Color', [0 0.3 0.8],    'LineWidth', 1.0);
    h3 = plot(ax, zz, S.frl, '-',  'Color', [0.85 0.2 0.1], 'LineWidth', 1.2);
    h4 = plot(ax, zz, S.fzl, '--', 'Color', [0.85 0.2 0.1], 'LineWidth', 1.0);
    xlim(ax, [-25 25]); xlabel(ax, 'z / mm'); ylabel(ax, 'f / (N·m^{-1}·A^{-1})');
    % 标签用英文
    legend(ax, [h1 h2 h3 h4], {'wall set, f_z', 'wall set, f_r (cross)', 'liquid set, f_r', 'liquid set, f_z (cross)'}, 'Location', 'southwest', 'Box', 'off', 'FontSize', 8);
    title(ax, sprintf('%s: wall set %d\\times%.2f mm (%d turns), liquid set %d\\times%.2f mm (%d turns)', S.id, numel(S.zw), ww_f*1e3, Nw_f, numel(S.zl), wl_f*1e3, Nl_f), ...
        'FontWeight', 'normal', 'FontSize', 9);
    set(ax, 'Box', 'on', 'Layer', 'top');
    ax = nexttile(tl, c + 2); hold(ax, 'on');
    ll = lam*1e3;
    h1 = plot(ax, ll, S.Fzw, '-',  'Color', [0 0.3 0.8],    'LineWidth', 1.2);
    h2 = plot(ax, ll, S.Frw, '--', 'Color', [0 0.3 0.8],    'LineWidth', 1.0);
    h3 = plot(ax, ll, S.Frl, '-',  'Color', [0.85 0.2 0.1], 'LineWidth', 1.2);
    h4 = plot(ax, ll, S.Fzl, '--', 'Color', [0.85 0.2 0.1], 'LineWidth', 1.0);
    yl = max([S.Fzw S.Frw S.Frl S.Fzl])*1.15;  ylim(ax, [0 yl]); xlim(ax, [2 30]);
    for j = 1:numel(lam_t)
        xline(ax, lam_t(j)*1e3, ':', sprintf('%.1f', lam_t(j)*1e3), 'Color', [0.4 0.4 0.4], ...
            'LabelOrientation', 'aligned', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'left', 'FontSize', 8);
    end
    xlabel(ax, '\lambda / mm'); ylabel(ax, '|F(\lambda)| / (N·A^{-1})');
    legend(ax, [h1 h2 h3 h4], {'wall set, F_z', 'wall set, F_r', 'liquid set, F_r', 'liquid set, F_z'}, 'Location', 'northeast', 'Box', 'off', 'FontSize', 8);
    title(ax, [S.id ' spatial spectra'], 'FontWeight', 'normal', 'FontSize', 9);
    set(ax, 'Box', 'on', 'Layer', 'top');
end
hideToolbars(fig2);
exportgraphics(fig2, fullfile(figdir, 'E2_2_final_spectra.png'), 'Resolution', 200);
savefig(fig2, fullfile(figdir, 'E2_2_final_spectra.fig'));

%% 8. 结果文本
fid = fopen(fullfile(here, 'E2_结果.txt'), 'w', 'n', 'UTF-8');
fprintf(fid, 'E2 线圈几何 结果  %s\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')));
fprintf(fid, '输入：E1_fields.mat（D0 / D1 / D1T 磁体叠列在 r_o = 6.5 mm 的场）；线圈导线中心半径 %.2f mm，0.15 mm 线单层；节距 10.4 mm。\n', geo.rc*1e3);
fprintf(fid, 'F(λ) 单位 N 每安培线圈电流（互易：也是接收灵敏度的度量）。目标 λ：10.4（壁 520 kHz）、6.3 / 5.2 / 4.5（液 300 kHz 三模态）。\n\n');
fprintf(fid, '== 壁组扫描（磁下元件，电流随磁体极性交替；列 = F_z(10.4)、F_r(10.4) 反串、F_z(5.2) 旁瓣、每匝 F_z(10.4)）==\n');
fprintf(fid, '%6s %4s', 'w_w', 'N');
for s = 1:3, fprintf(fid, ' | %7s %7s %7s %8s', [St{s}.id ' Fz'], 'Fr_x', 'Fz5.2', 'Fz/N'); end
fprintf(fid, '\n');
for i = 1:numel(ww)
    fprintf(fid, '%6.2f %4d', ww(i)*1e3, St{1}.Nw(i));
    for s = 1:3, fprintf(fid, ' | %7.3f %7.3f %7.3f %8.4f', St{s}.Fw(i,1), St{s}.Fw_x(i,1), St{s}.Fw(i,3), St{s}.Fw(i,1)/St{s}.Nw(i)); end
    fprintf(fid, '\n');
end
fprintf(fid, '\n== 液组扫描（缝下元件，两缝电流反向；列 = F_r(6.3)、F_r(5.2)、F_r(4.5)、F_z(5.2) 反串、每匝 F_r(5.2)）==\n');
fprintf(fid, '%6s %4s', 'w_l', 'N');
for s = 1:3, fprintf(fid, ' | %7s %7s %7s %7s %8s', [St{s}.id ' 6.3'], '5.2', '4.5', 'Fz_x', 'Fr/N'); end
fprintf(fid, '\n');
for i = 1:numel(wl)
    fprintf(fid, '%6.2f %4d', wl(i)*1e3, St{1}.Nl(i));
    for s = 1:3, fprintf(fid, ' | %7.3f %7.3f %7.3f %7.3f %8.4f', St{s}.Fl(i,2), St{s}.Fl(i,3), St{s}.Fl(i,4), St{s}.Fl_x(i,3), St{s}.Fl(i,3)/St{s}.Nl(i)); end
    fprintf(fid, '\n');
end
fprintf(fid, '\n== 各叠列自己的最优宽（壁组按 F_z(10.4) 最大，液组按 F_r(5.2) 最大）==\n');
for s = 1:3
    fprintf(fid, '%-4s 壁组 w_w = %.2f mm（%d 匝）F_z(10.4) = %.3f；液组 w_l = %.2f mm（%d 匝）F_r(5.2) = %.3f\n', ...
        St{s}.id, St{s}.ww_opt*1e3, St{s}.Nw_opt, St{s}.Fw(St{s}.iw,1), St{s}.wl_opt*1e3, St{s}.Nl_opt, St{s}.Fl(St{s}.il,3));
end
fprintf(fid, '\n== 统一规格：壁组元件 %.2f mm、%d 匝单层（λ_w/2；F_z(10.4) 为各叠列最大值的 %.0f%% / %.0f%% / %.0f%%）；液组元件 %.2f mm、%d 匝单层；节距 10.4 mm ==\n', ...
    ww_f*1e3, Nw_f, 100*ratio_w(1), 100*ratio_w(2), 100*ratio_w(3), wl_f*1e3, Nl_f);
fprintf(fid, '%-4s | %8s %8s %6s | %8s %8s %8s %8s %6s | %10s %10s\n', 'id', 'Fzw10.4', 'Frw10.4', 'ax dB', 'Frl6.3', 'Frl5.2', 'Frl4.5', 'Fzl5.2', 'rd dB', '单套Fz5.2', '单套Fr5.2');
for s = 1:3
    S = St{s};
    fprintf(fid, '%-4s | %8.3f %8.3f %6.1f | %8.3f %8.3f %8.3f %8.3f %6.1f | %10.3f %10.3f\n', S.id, ...
        S.Fzw_t(1), S.Frw_t(1), 20*log10(S.Fzw_t(1)/S.Frw_t(1)), S.Frl_t(2), S.Frl_t(3), S.Frl_t(4), S.Fzl_t(3), ...
        20*log10(S.Frl_t(3)/S.Fzl_t(3)), S.Fz_single_t(3), S.Fr_single_t(3));
end
fprintf(fid, '"单套"列：磁下 2.6 mm（17 匝）窄元件在 λ = 5.2 的轴向力 F_z 与径向力 F_r——方向仍以轴向为主，不能替代缝下的液组。\n');
fprintf(fid, '壁组径向反串的作用：叠列对称（f_z 偶、f_r 奇），两个传播方向幅值相同，无方向性；S0 外表面 u_r 由泊松变薄给出、与 u_z 正交相位，\n');
fprintf(fid, '  反串项与轴向项同号相加。薄板近似 η = ν h k/(2(1−ν)) = %.3f（真实 S0 表面比约 0.2–0.3），增益 ≈ 1 + η F_r/F_z：D1T 壁组 %.2f–%.2f（E3 核）。\n', ...
    eta_plate, 1 + eta_plate*St{3}.Frw_t(1)/St{3}.Fzw_t(1), 1 + 0.3*St{3}.Frw_t(1)/St{3}.Fzw_t(1));
fprintf(fid, '\n== 电感与电阻粗估（Wheeler 空心；管壁镜像再降 30–50%%，以 LCR 表实测为准）==\n');
fprintf(fid, '壁组单元件 %d 匝：L ≈ %.1f µH 空心（估 %.1f–%.1f µH 装上管），R_dc = %.2f Ω；%d 元件串联 ×%d\n', ...
    Nw_f, L_w1, 0.5*L_w1, 0.7*L_w1, Rdc(Nw_f), 3, 3);
fprintf(fid, '液组单元件 %d 匝：L ≈ %.2f µH 空心（估 %.2f–%.2f µH 装上管），R_dc = %.2f Ω；2 元件串联 ×2\n', ...
    Nl_f, L_l1, 0.5*L_l1, 0.7*L_l1, Rdc(Nl_f));
fprintf(fid, '双层绕制：匝数 ×2、F ×2、L ×4；壁组 3 元件双层 |Z|(520 kHz) ≈ %.0f Ω，液组 2 元件双层 |Z|(300 kHz) ≈ %.0f Ω（装上管估值取 0.6 倍空心）。\n', ...
    2*pi*520e3*0.6*4*3*L_w1*1e-6, 2*pi*300e3*0.6*4*2*L_l1*1e-6);
fclose(fid);

%% 9. 给 E3/E4 的线圈规格
coils.pitch_mm = 10.4;  coils.wire_mm = 0.15;  coils.rc_mm = geo.rc*1e3;
coils.wall   = struct('width_mm', ww_f*1e3, 'turns', Nw_f, 'z_mm_D0', St{1}.zw*1e3, 'sign_D0', St{1}.sw, ...
                      'z_mm_D1', St{2}.zw*1e3, 'sign_D1', St{2}.sw, 'L_uH_air', L_w1);
coils.liquid = struct('width_mm', wl_f*1e3, 'turns', Nl_f, 'z_mm_D0', St{1}.zl*1e3, 'sign_D0', St{1}.sl, ...
                      'z_mm_D1', St{2}.zl*1e3, 'sign_D1', St{2}.sl, 'L_uH_air', L_l1);
coils.sweep = struct('ww', ww, 'wl', wl, 'St', {St});
save(fullfile(here, 'E2_coils.mat'), '-struct', 'coils');

figure(fig1); figure(fig2); drawnow;
fprintf('E2 完成：壁组 %.2f mm / %d 匝，液组 %.2f mm / %d 匝\n', ww_f*1e3, Nw_f, wl_f*1e3, Nl_f);
end

%% ---------- 局部函数 ----------
function Fv = spec(f, z, lam)
k  = 2*pi./lam(:).';
Fv = abs(trapz(z, f.*exp(-1i*z*k), 1));
end

function hideToolbars(fig)
axs = findall(fig, 'Type', 'axes');
for k = 1:numel(axs)
    try, axs(k).Toolbar.Visible = 'off'; catch, end
end
end
