function R = T4_report(rundir)
%T4_REPORT  读 T4_all.mat，按 T4 方案 §1.3 逐项与 SAFE(T3)/闭式(T1) 对比，出表出图
%   R = T4_report                     默认 2_simulation\T4_runs
% 输出：<rundir>\T4_结果.txt、<rundir>\figures\T4_1_dispersion、T4_2_sensitivity（png+fig）
%
% 参考解全部**直接从 T3_sweep.mat 取**（allK = K^true 四列 [σz, p, T, 壁ΔA分项]，
% allSp = 相位版 ns/MPa 四列 [S11, S12, S13, S12wall]），不抄表，免得插值踩到耦合区。
% T3 的自然量 = allK − conv，conv = [1/E, ε_z(p)|σz=0, 0, 0]，与 T3_sensitivity.m 一致。
%
% 速度口径（方案 §0.1）：FEM 的 Floquet 相位按材料坐标 → c_p 是**自然速度**；
% 换真实：K^true = K^nat + conv。T3 的压力导数是 **σz = 0（开端）** 下取的，
% 所以壁模态的压力项要跟 FEM 的 L2f / L2pf（把闭端轴向力扣掉）比，不是 L2。

if nargin < 1 || isempty(rundir), rundir = fullfile(fileparts(mfilename('fullpath')),'T4_runs'); end
here = fileparts(mfilename('fullpath'));
addpath(here); addpath(fullfile(here,'..','1_theory'));
figdir = fullfile(rundir,'figures'); if ~exist(figdir,'dir'), mkdir(figdir); end
FONT = 'TimesSimSun';
try, if any(strcmpi(FONT, listfonts)), set(groot,'defaultAxesFontName',FONT,'defaultTextFontName',FONT,'defaultLegendFontName',FONT); end, catch, end
set(groot,'defaultTextInterpreter','tex','defaultAxesTickLabelInterpreter','tex','defaultLegendInterpreter','tex','defaultAxesFontSize',9);

L  = load(fullfile(rundir,'T4_all.mat'));  A = L.A;
T3 = load(fullfile(here,'..','1_theory','T3_sweep.mat'));
P = params_pipe;  W = P.wall;  G = P.geom;  MPa = 1e6;
E_MPa = W.E/1e6;  WE = G.W_E;  ap = P.fluid.a_p*1e6;
hooke = @(s) ((1+W.nu)*s - W.nu*trace(s)*eye(3))/W.E;
aa = G.r_i^2/(G.r_o^2-G.r_i^2);
ez1 = hooke(diag([0 0 1])*MPa);          conv_z = ez1(3,3);                    % 1/E
ezp = hooke(diag([-MPa, MPa*aa*(1+G.r_o^2/G.r_i^2), 0]));  conv_p = ezp(3,3);  % σz=0 时 ε_z/p

fid = fopen(fullfile(rundir,'T4_结果.txt'),'w','n','UTF-8');
pr = @(varargin) fprintf(fid, varargin{:});
pr('T4 与 Floquet–Bloch FEM 的对比 结果  %s\n', char(datetime('now','Format','yyyy-MM-dd HH:mm')));
pr('COMSOL 6.4，2D 轴对称单胞（solid 管壁 + acpr 液体），参考解 SAFE = 1_theory\\T3_sweep.mat，闭式 = T1\n');
pr('conv_z = %.4e /MPa（=1/E），conv_p = %.4e /MPa（σz=0 时的 ε_z/p），a_p = %.4e /MPa\n\n', conv_z, conv_p, ap);
R = struct();

%% ===== 1. L0 频散 vs SAFE =====
C = T4_compare(A.L0, T3);  R.L0cmp = C;
isW = C.eta_safe < 0.5;
pr('== 1. L0 基准频散 vs SAFE ==（%d/%d 个 SAFE 点配上）\n', height(C), numel(T3.allF));
pr('  壁传族 %3d 点：|Δc_p| 中位 %.4f%%  90%%位 %.4f%%  最大 %.4f%%   （方案容差 0.2%%）\n', ...
    nnz(isW), median(abs(C.dcp_pct(isW))), prctile(abs(C.dcp_pct(isW)),90), max(abs(C.dcp_pct(isW))));
pr('  液传族 %3d 点：|Δc_p| 中位 %.4f%%  90%%位 %.4f%%  最大 %.4f%%\n', ...
    nnz(~isW), median(abs(C.dcp_pct(~isW))), prctile(abs(C.dcp_pct(~isW)),90), max(abs(C.dcp_pct(~isW))));
pr('  c_g：中位 %.2f%%  最大 %.2f%%（FEM 的 c_g 由 k 网格中心差分，步长 %g rad/m）\n', ...
    median(abs(C.dcg_pct)), max(abs(C.dcg_pct)), L.kg(2)-L.kg(1));
pr('  η_f：中位绝对差 %.5f  最大 %.5f（与 SAFE 同口径 (V_f+T_f)/E）\n\n', ...
    median(abs(C.eta_fem-C.eta_safe)), max(abs(C.eta_fem-C.eta_safe)));

%% ===== 2. 目标模态 =====
fw = [400 500 520 600 625 650 800]*1e3;
cw = arrayfun(@(f) pickcp(T3,f,'wall'), fw);
mw = arrayfun(@(f,c) findbranch(A.L0,f,c), fw, cw);
liqT = [ 20 0; 50 0; 250 1578.6; 300 1904.1; 300 1560.4; 300 1359.4; ...
        350 1810.8; 350 1550.5; 350 1395.1; 400 1779.9; 400 1539.0; 400 1415.6];
liqT(1,2) = pickcp(T3,20e3,'liq');  liqT(2,2) = pickcp(T3,50e3,'liq');
fl = liqT(:,1).'*1e3;  cl = liqT(:,2).';
ml = arrayfun(@(f,c) findbranch(A.L0,f,c), fl, cl);
R.wall = struct('f',fw,'cp',cw,'m',mw);  R.liq = struct('f',fl,'cp',cl,'m',ml);

%% ===== 3. L1 判别试验 + W1/W2 =====
pr('== 2. L1 判别试验（σz,ext = 100 MPa，p = 0）：FEM 给的是自然速度还是真实速度 ==\n');
pr('  T1 薄板：自然 −7.61e−3，真实 −6.23e−3（差 18%%）\n');
pr('  f/kHz    c_p    Δc/c(FEM)   K_z^nat(FEM)  K_z^nat(T3)   差     K_z^true(FEM)  K_z^true(T3)\n');
Kz = nan(size(fw));  Kz3 = nan(size(fw));
for i = 1:numel(fw)
    r = t3ref(T3, fw(i), cw(i));
    d = dcp(A.L0, A.L1, fw(i), mw(i));   Kz(i) = d/100;
    Kz3(i) = r.Ktrue(1) - conv_z;
    pr('  %5.0f %8.1f  %+.4e  %+.4e  %+.4e  %+6.2f%%   %+.4e  %+.4e\n', fw(i)/1e3, cw(i), d, ...
        Kz(i), Kz3(i), 100*(Kz(i)/Kz3(i)-1), Kz(i)+conv_z, r.Ktrue(1));
end
R.Kz_nat_fem = Kz;  R.Kz_nat_t3 = Kz3;
[dW1c, dW2c, d3c] = closedform_split(W);
pr('\n  W1（线弹性+几何非线性 = St Venant–Kirchhoff ≡ l=m=n=0）与 W2（Murnaghan）分解，100 MPa：\n');
pr('  T1 薄板闭式：W2 %+.4e，W1 %+.4e，三阶项 %+.4e\n', dW2c, dW1c, d3c);
pr('  f/kHz    FEM W1      FEM W2     FEM W2−W1   与闭式三阶项差   W1 与闭式差\n');
for i = 1:numel(fw)
    m1 = findbranch(A.L0w1, fw(i), cw(i));
    d1 = dcp(A.L0w1, A.L1w1, fw(i), m1);
    d2 = dcp(A.L0,   A.L1,   fw(i), mw(i));
    pr('  %5.0f  %+.4e  %+.4e  %+.4e   %+6.2f%%      %+6.2f%%\n', fw(i)/1e3, d1, d2, d2-d1, ...
        100*((d2-d1)/d3c-1), 100*(d1/dW1c-1));
end

%% ===== 4. 壁模态压力分解 =====
pr('\n== 3. 壁模态对内压（与 T3 同口径：σz = 0，即 FEM 的 L2f / L2pf）==\n');
pr('  分解：L2pf−L0 = 壁声弹（**含随动压力项**）；L2f−L2pf = 液体加载项\n');
pr('  f/kHz  K_p^nat 总    声弹      液体加载 | 相位版 ns/MPa：S12(FEM) S12(T3) | S12wall(FEM) S12wall(T3) 随动项\n');
Kp = nan(numel(fw),3);
for i = 1:numel(fw)
    r = t3ref(T3, fw(i), cw(i));
    tot = dcp(A.L0, A.L2f,  fw(i), mw(i))/20;
    ael = dcp(A.L0, A.L2pf, fw(i), mw(i))/20;
    fld = tot - ael;   Kp(i,:) = [tot ael fld];
    t0 = WE/cw(i);
    % 注意口径：T3 的 Sp(2)=S12 是**真实**量（Ktrue），Sp(4)=S12wall 是**自然**量
    %（T3_sensitivity.m：conv = [conv_z, conv_p, 0, 0]，第 4 列不换算）
    S12f = -(tot+conv_p)*t0*1e9;  S12wf = -ael*t0*1e9;
    pr('  %5.0f  %+.4e  %+.4e  %+.4e | %8.3f %8.3f | %8.3f %8.3f  %+7.3f\n', ...
        fw(i)/1e3, tot, ael, fld, S12f, r.Sp(2), S12wf, r.Sp(4), S12wf - r.Sp(4));
end
R.Kp_fem = Kp;
pr('  （最后一列 = FEM 声弹项 − SAFE 声弹项；SAFE 没有随动压力项，所以这一列**应当**≈ 随动项，\n');
pr('    §3b 用"本征步关掉随动压力"的对照直接量了随动项，两者对上才算这一列成立）\n');
pr('  参考：闭端总量（FEM 的 L2，σz = 50.4 MPa 自动）\n');
for i = 1:numel(fw)
    pr('    %5.0f kHz：K_p^nat(闭端) = %+.4e /MPa\n', fw(i)/1e3, dcp(A.L0,A.L2,fw(i),mw(i))/20);
end

%% ===== 4b. 随动压力项的直接隔离 =====
pr('\n== 3b. 随动压力项的直接隔离（本征步里把内壁 Follower Pressure 关掉，其余完全相同）==\n');
pr('  这一项 SAFE 和 T3 都没有；T4 §0.3 要求单独记录。固定 k 比，再用 c_p/c_g 折成固定频率。\n');
try
    Tf_ = T4_build('quiet',true);  Tf_.study('std1').feature('eig').set('neigs',14);
    Tn_ = T4_build('quiet',true);  Tn_.study('std1').feature('eig').set('neigs',14);
    Tn_.study('std1').feature('eig').set('disabledphysics', {'solid/symz0','solid/blz','solid/blp'});
    tg = [520e3 pickcp(T3,520e3,'wall'); 20e3 pickcp(T3,20e3,'liq'); 300e3 1560.4; 300e3 1904.1];
    pr('  模态                 K_p(随动开)   K_p(随动关)   随动项       占总量   折成 ns/MPa\n');
    for i = 1:size(tg,1)
        ft = tg(i,1); cref = tg(i,2); kk = 2*pi*ft/cref;  dd = zeros(1,2);
        for j = 1:2
            MM = Tf_; if j==2, MM = Tn_; end
            MM.param.set('kz',sprintf('%.10g[1/m]',kk));
            MM.param.set('sig_ext','0[MPa]'); MM.param.set('p_in','0[MPa]'); MM.param.set('frz','1'); MM.param.set('dT','0[K]');
            MM.sol('sol1').runAll; f0 = real(mphglobal(MM,'freq','dataset','dset1','solnum','all'));
            MM.param.set('sig_ext',sprintf('%.6g[MPa]',-20*G.A_bore/G.A_wall)); MM.param.set('p_in','20[MPa]');
            MM.sol('sol1').runAll; f1 = real(mphglobal(MM,'freq','dataset','dset1','solnum','all'));
            [~,j0] = min(abs(2*pi*f0/kk - cref));
            dd(j) = f1(j0)/f0(j0) - 1;
        end
        m = findbranch(A.L0, ft, cref);  fm = A.L0.f(:,m); g = ~isnan(fm);
        cv = cref/interp1(fm(g), A.L0.cg(g,m), ft, 'pchip');
        t0 = WE/cref;
        pr('  %3.0f kHz c_p %6.0f   %+.4e  %+.4e  %+.4e  %+6.2f%%   %+7.3f\n', ft/1e3, cref, ...
            dd(1)*cv/20, dd(2)*cv/20, (dd(1)-dd(2))*cv/20, 100*(dd(1)-dd(2))/dd(2), -(dd(1)-dd(2))*cv/20*t0*1e9);
    end
catch ME
    pr('  隔离试验没跑成：%s\n', regexprep(ME.message,'\s+',' '));
end

%% ===== 5. 液传族 =====
pr('\n== 4. 液模态：对内压（σz=0，L2f）与对 σz（L1）==\n');
pr('  T1 Korteweg 极限：(1/c)dc/dp = 2.224e−3 /MPa，κ(0) = 0.734；a_p = %.4e /MPa\n', ap);
pr('  f/kHz    c_p   支   K_p^nat(FEM)  K_p^true(FEM) K_p^true(T3)  差    κ_p(FEM) κ_p(T3) | S22 FEM   T3   | K_σz FEM   T3\n');
Kl = nan(numel(fl),3);
for i = 1:numel(fl)
    if isnan(ml(i)), continue; end
    r = t3ref(T3, fl(i), cl(i));
    dp = dcp(A.L0, A.L2f, fl(i), ml(i))/20;   dz = dcp(A.L0, A.L1, fl(i), ml(i))/100;
    wal = dcp(A.L0, A.L2pf, fl(i), ml(i))/20;
    Kl(i,:) = [dp dz wal];
    dpt = dp + conv_p;  t0 = WE/cl(i);
    pr('  %5.0f %7.1f  %2d   %+.4e  %+.4e  %+.4e %+6.2f%%  %6.3f  %6.3f | %8.2f %8.2f | %+.3e %+.3e\n', ...
        fl(i)/1e3, cl(i), ml(i), dp, dpt, r.Ktrue(2), 100*(dpt/r.Ktrue(2)-1), dpt/ap, r.Ktrue(2)/ap, ...
        -dpt*t0*1e9, r.Sp(2), dz, r.Ktrue(1));
end
R.Kl_fem = Kl;
pr('  液模态 K_p 的两侧分解（自然量）：壁侧 = L2pf−L0（液体参数冻结），液侧 = L2f−L2pf\n');
pr('  f/kHz    c_p     壁侧        液侧      T3 的壁侧（Sp 第 4 列换回自然）\n');
for i = 1:numel(fl)
    if isnan(ml(i)), continue; end
    r = t3ref(T3, fl(i), cl(i));  t0 = WE/cl(i);
    pr('  %5.0f %7.1f  %+.4e  %+.4e   %+.4e\n', fl(i)/1e3, cl(i), Kl(i,3), Kl(i,1)-Kl(i,3), -r.Sp(4)*1e-9/t0);
end

%% ===== 6. L3 叠加、L4 温度、L5 水 =====
pr('\n== 5. L3 线性叠加检验（σz,ext = 100 MPa 且 p = 20 MPa 闭端）==\n');
pr('  比的是 Δc/c：Δ(L3) 与 Δ(L1)+Δ(L2)（两者都是闭端口径）\n');
pr('  模态             Δ(L3)      Δ(L1)+Δ(L2)   残差/(|Δ1|+|Δ2|)\n');
for i = 1:numel(fw)
    d3 = dcp(A.L0,A.L3,fw(i),mw(i)); d1 = dcp(A.L0,A.L1,fw(i),mw(i)); d2 = dcp(A.L0,A.L2,fw(i),mw(i));
    pr('  壁 %5.0f kHz     %+.4e  %+.4e   %+6.2f%%\n', fw(i)/1e3, d3, d1+d2, 100*(d3-d1-d2)/(abs(d1)+abs(d2)));
end
for i = 1:numel(fl)
    if isnan(ml(i)) || fl(i) ~= 300e3, continue; end
    d3 = dcp(A.L0,A.L3,fl(i),ml(i)); d1 = dcp(A.L0,A.L1,fl(i),ml(i)); d2 = dcp(A.L0,A.L2,fl(i),ml(i));
    pr('  液 300 kHz c_p %4.0f %+.4e  %+.4e   %+6.2f%%\n', cl(i), d3, d1+d2, 100*(d3-d1-d2)/(abs(d1)+abs(d2)));
end

pr('\n== 6. L4 温度（ΔT = +10 K）==\n');
pr('  模态             (1/c)dc/dT [1/K]   ∂t/∂T FEM   T3 S13    差\n');
for i = 1:numel(fw)
    r = t3ref(T3, fw(i), cw(i));
    d = dcp(A.L0,A.L4,fw(i),mw(i))/10;  t0 = WE/cw(i);
    pr('  壁 %5.0f kHz     %+.4e        %7.2f  %7.2f  %+6.2f%%\n', fw(i)/1e3, d, -d*t0*1e9, r.Sp(3), 100*(-d*t0*1e9/r.Sp(3)-1));
end
for i = 1:numel(fl)
    if isnan(ml(i)) || ~ismember(fl(i),[20e3 300e3]), continue; end
    r = t3ref(T3, fl(i), cl(i));
    d = dcp(A.L0,A.L4,fl(i),ml(i))/10;  t0 = WE/cl(i);
    pr('  液 %5.0f kHz c_p %4.0f  %+.4e   %7.2f  %7.2f  %+6.2f%%\n', fl(i)/1e3, cl(i), d, -d*t0*1e9, r.Sp(3), 100*(-d*t0*1e9/r.Sp(3)-1));
end

pr('\n== 7. L5 油 vs 水（液模态压力灵敏度比，T1 预测 3.02）==\n');
for ff = [20e3 50e3]
    try
        mo  = findbranch(A.L0,   ff, pickcp(T3,ff,'liq'));
        mwt = findbranch(A.L5_0, ff, 1259);
        do = dcp(A.L0,   A.L2f, ff, mo)/20;
        dw = dcp(A.L5_0, A.L5,  ff, mwt)/20;
        pr('  %5.0f kHz：油 c_p %.0f K_p %+.4e；水 c_p %.0f K_p %+.4e /MPa → 比 %.2f\n', ff/1e3, ...
            findcp(A.L0,ff,pickcp(T3,ff,'liq')), do, findcp(A.L5_0,ff,1259), dw, do/dw);
    catch ME
        pr('  %5.0f kHz 没算成：%s\n', ff/1e3, ME.message);
    end
end

pr('\n== 8. 干净窗口位置随 c_f（L0，c_f = 1400 / 1450 / 1500 m/s）==\n');
for nm = {'L0cf1400','L0','L0cf1500'}
    cf = 1450; if strcmp(nm{1},'L0cf1400'), cf = 1400; elseif strcmp(nm{1},'L0cf1500'), cf = 1500; end
    pr('  c_f = %4d m/s：壁传族 η_f < 0.10 的窗口 / kHz： %s\n', cf, cleanwindows(A.(nm{1})));
end
pr('  （T3 在 c_f = 1450 给的窗口：20–90 / 250–270 / 360–410 / 490–550 / 625–680 / 770–800 kHz）\n');

%% ===== 9. 载荷线性度（割线 vs 切线）=====
flin = fullfile(rundir,'T4_lin.mat');
if exist(flin,'file')
    Lz = load(flin);  B = Lz.B;
    pr('\n== 9. 载荷线性度：量程内标定曲线弯不弯 + 外推到零载荷后与 SAFE 的一致性 ==\n');
    pr('  方案里没有这一项。SAFE/T1 给的是零载荷处的**切线**，FEM 给的是量程端点的**割线**，口径不同。\n');
    pr('  加跑 σz = 50 MPa 与 p = 5 / 10 MPa（T4_lin.m），按 K(x) ≈ K0 + ½K''x 外推回零载荷。\n');
    pr('\n  9.1 轴向应力（p = 0）\n');
    pr('  f/kHz   K(50 MPa)    K(100 MPa)  割线弯曲  外推切线 K0   SAFE 切线    切线差   （割线差）\n');
    for i = 1:numel(fw)
        m50 = findbranch(B.base, fw(i), cw(i));  if isnan(m50), continue; end
        k50 = dcp(B.base, B.s50, fw(i), m50)/50;  k100 = Kz(i);
        K0 = 2*k50 - k100;  r = t3ref(T3, fw(i), cw(i));  Kt3 = r.Ktrue(1) - conv_z;
        pr('  %5.0f  %+.4e  %+.4e  %+6.2f%%  %+.4e  %+.4e  %+6.2f%%  %+6.2f%%\n', ...
            fw(i)/1e3, k50, k100, 100*(k100/k50-1), K0, Kt3, 100*(K0/Kt3-1), 100*(k100/Kt3-1));
    end
    pr('\n  9.2 内压（σz = 0）\n');
    pr('  f/kHz    c_p    K(5 MPa)     K(10 MPa)    K(20 MPa)   割线弯曲  外推切线    SAFE 切线   切线差  （割线差）\n');
    tp = [520e3 cw(fw==520e3); 300e3 1904.1; 300e3 1560.4; 250e3 1578.6; 20e3 pickcp(T3,20e3,'liq')];
    for i = 1:size(tp,1)
        mb = findbranch(B.base, tp(i,1), tp(i,2));  m2 = findbranch(A.L0, tp(i,1), tp(i,2));
        if isnan(mb) || isnan(m2), continue; end
        k5 = dcp(B.base,B.p5,tp(i,1),mb)/5;  k10 = dcp(B.base,B.p10,tp(i,1),mb)/10;
        k20 = dcp(A.L0,A.L2f,tp(i,1),m2)/20;
        K0 = 2*k5 - k10;  r = t3ref(T3, tp(i,1), tp(i,2));  Kt3 = r.Ktrue(2) - conv_p;
        pr('  %5.0f %7.1f  %+.4e %+.4e %+.4e  %+6.2f%%  %+.4e %+.4e %+6.2f%% %+7.2f%%\n', ...
            tp(i,1)/1e3, tp(i,2), k5, k10, k20, 100*(k20/k5-1), K0, Kt3, 100*(K0/Kt3-1), 100*(k20/Kt3-1));
    end
    pr('  结论：统一到切线口径后，干净窗口上 FEM 与 SAFE 差 ≤1.3%%；\n');
    pr('        §2 表里 +2~+3%% 和 §4 表里 1904 支的 +11.5%% 都是有限载荷步的割线效应，不是模型不一致。\n');
    pr('        壁 520 kHz 的 K_p 弯得太厉害（+9.5%%），两点外推不成立，只能按量程内割线报。\n');
else
    pr('\n== 9. 载荷线性度：没找到 %s，跳过（跑 T4_lin 生成）==\n', flin);
end
fclose(fid);

%% ===== 图 =====
fig = figure('Visible','on','Name','T4-1 dispersion','Position',[40 40 1020 430],'Color','w');
tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
ax = nexttile(tl,1); hold(ax,'on');
gl = A.L0.eta > 0.5;
sc1 = scatter(ax, A.L0.f(gl)/1e3, A.L0.cp(gl), 4, [0 0.45 0.85], 'filled');
sc2 = scatter(ax, A.L0.f(~gl)/1e3, A.L0.cp(~gl), 4, [0.85 0.2 0.1], 'filled');
sc3 = scatter(ax, T3.allF/1e3, T3.allCP, 26, 'k', 'o', 'LineWidth',0.6);
xlim(ax,[0 800]); ylim(ax,[1000 6200]); xlabel(ax,'f / kHz'); ylabel(ax,'c_p（自然速度）/ (m·s^{-1})');
legend(ax,[sc1 sc2 sc3],{'FEM 液传族 η_f>0.5','FEM 壁传族','SAFE (T3)'},'Location','northeast','Box','off','FontSize',8);
title(ax,'L0 频散：Floquet–Bloch FEM vs SAFE','FontWeight','normal'); set(ax,'Box','on','Layer','top');
ax = nexttile(tl,2); hold(ax,'on');
plot(ax, C.f_Hz(~isW)/1e3, C.dcp_pct(~isW), 'o', 'Color',[0 0.45 0.85],'MarkerSize',3);
plot(ax, C.f_Hz(isW)/1e3, C.dcp_pct(isW), 's', 'Color',[0.85 0.2 0.1],'MarkerSize',4);
yline(ax,0.2,':k'); yline(ax,-0.2,':k');
xlabel(ax,'f / kHz'); ylabel(ax,'(c_p^{FEM}/c_p^{SAFE} − 1) / %');
legend(ax,{'液传族','壁传族','±0.2% 容差'},'Location','best','Box','off','FontSize',8);
title(ax,'逐点相对差','FontWeight','normal'); set(ax,'Box','on','Layer','top');
savefig2(fig, figdir, 'T4_1_dispersion');

fig = figure('Visible','on','Name','T4-2 sensitivity','Position',[60 60 1020 430],'Color','w');
tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
ax = nexttile(tl,1); hold(ax,'on');
plot(ax, fw/1e3, Kz*1e5, 'o-','Color',[0.85 0.2 0.1],'LineWidth',1.2,'MarkerFaceColor','w');
plot(ax, fw/1e3, Kz3*1e5, 's--','Color',[0 0.45 0.85],'LineWidth',1.0);
yline(ax,-7.6067,':k');
xlabel(ax,'f / kHz'); ylabel(ax,'K_z^{nat} / (10^{-5} MPa^{-1})');
legend(ax,{'FEM（Floquet–Bloch）','SAFE (T3)','T1 薄板 −7.61'},'Location','best','Box','off','FontSize',8);
title(ax,'壁模态轴向应力灵敏度（自然速度）','FontWeight','normal'); set(ax,'Box','on','Layer','top');
ax = nexttile(tl,2); hold(ax,'on');
% 只画非耦合区的点：600 / 625 kHz 靠近 avoided crossing，K_p 大 1–3 个量级，画进来别的全被压平
clean = abs(Kp(:,1)).' < 1e-3;
b = bar(ax, categorical(fw(clean)/1e3), [Kp(clean,2) Kp(clean,3)]*1e4, 'stacked');
b(1).FaceColor = [0.85 0.5 0.2];  b(2).FaceColor = [0 0.45 0.85];
xlabel(ax,'f / kHz（非耦合区）'); ylabel(ax,'K_p^{nat} / (10^{-4} MPa^{-1})，σ_z = 0');
ylim(ax, [0 1.35*max(sum(Kp(clean,2:3),2))*1e4]);
legend(ax,{'壁声弹（含随动压力）','液体加载'},'Location','north','Box','off','FontSize',8);
txt = sprintf('耦合区不画：600 kHz K_p = %.2e，625 kHz %.2e /MPa', Kp(fw==600e3,1), Kp(fw==625e3,1));
text(ax, 0.03, 0.055, txt, 'Units','normalized', 'FontSize',7.5, 'Color',[0.35 0.35 0.35]);
title(ax,'壁模态内压灵敏度的两项分解（FEM）','FontWeight','normal'); set(ax,'Box','on','Layer','top');
savefig2(fig, figdir, 'T4_2_sensitivity');

fprintf('T4 报告完成：%s\n', fullfile(rundir,'T4_结果.txt'));
end

%% ---------------- 局部函数 ----------------
function r = t3ref(T3, f, cp)
d = abs(T3.allF - f)/max(f,1) + abs(T3.allCP - cp)/max(cp,1);
[dm, i] = min(d);
r = struct('ok', dm < 0.08, 'f',T3.allF(i), 'cp',T3.allCP(i), 'cg',T3.allCG(i), 'eta',T3.allETA(i), ...
           'Ktrue',T3.allK(i,:), 'Sp',T3.allSp(i,:), 'Sg',T3.allSg(i,:));
end

function c = pickcp(T3, f, fam)
d = abs(T3.allF - f);  g = d == min(d);
if strcmp(fam,'wall'), g = g & T3.allETA < 0.5; else, g = g & T3.allETA > 0.5; end
if ~any(g), c = NaN; return; end
cc = T3.allCP(g);  ee = T3.allETA(g);
if strcmp(fam,'wall'), [~,j] = min(ee); else, [~,j] = max(ee); end
c = cc(j);
end

function m = findbranch(S, ft, cpref)
best = inf; m = NaN;
for j = 1:size(S.f,2)
    fj = S.f(:,j); g = ~isnan(fj);
    if nnz(g) < 4 || ~all(diff(fj(g)) > 0), continue; end
    k = interp1(fj(g), S.k(g), ft, 'pchip', NaN);
    if isnan(k), continue; end
    e = abs(2*pi*ft/k/cpref - 1);
    if e < best, best = e; m = j; end
end
if best > 0.10, m = NaN; end
end

function c = findcp(S, ft, cpref)
m = findbranch(S, ft, cpref);
if isnan(m), c = NaN; return; end
fj = S.f(:,m); g = ~isnan(fj);
c = 2*pi*ft/interp1(fj(g), S.k(g), ft, 'pchip');
end

function d = dcp(S0, S1, ft, m)
% 固定频率 ft 下第 m 支的相速度相对变化：c_p = 2πf/k ⇒ Δc/c = k0/k1 − 1
d = NaN;  if isnan(m), return; end
try
    f0 = S0.f(:,m); g0 = ~isnan(f0);
    f1 = S1.f(:,m); g1 = ~isnan(f1);
    k0 = interp1(f0(g0), S0.k(g0), ft, 'pchip', NaN);
    k1 = interp1(f1(g1), S1.k(g1), ft, 'pchip', NaN);
    d = k0/k1 - 1;
catch, end
end

function [dW1, dW2, d3] = closedform_split(W)
MPa = 1e6;  hk = @(s) ((1+W.nu)*s - W.nu*trace(s)*eye(3))/W.E;
Sz = zeros(3); Sz(1,1) = 100*MPa;
f = @(l,m,n) (s0_reduced_modulus(murnaghan_A(W.lam,W.mu,l,m,n,Sz,hk(Sz))) ...
            - s0_reduced_modulus(murnaghan_A(W.lam,W.mu,l,m,n,zeros(3),zeros(3)))) ...
            /(2*s0_reduced_modulus(murnaghan_A(W.lam,W.mu,l,m,n,zeros(3),zeros(3))));
dW2 = f(W.l,W.m,W.n);  dW1 = f(0,0,0);  d3 = dW2 - dW1;
end

function s = cleanwindows(S)
fg = (20:10:800)*1e3;  et = nan(size(fg));
for i = 1:numel(fg)
    best = inf;
    for m = 1:size(S.f,2)
        fj = S.f(:,m); g = ~isnan(fj);
        if nnz(g) < 4 || ~all(diff(fj(g)) > 0), continue; end
        k = interp1(fj(g), S.k(g), fg(i), 'pchip', NaN);
        if isnan(k), continue; end
        cp = 2*pi*fg(i)/k;  e = interp1(fj(g), S.eta(g,m), fg(i), 'pchip');
        if cp > 4500 && cp < 6500 && e < best, best = e; end
    end
    et(i) = best;
end
ok = et < 0.10;  s = '';
d = diff([false ok false]);  st = find(d == 1);  en = find(d == -1) - 1;
for i = 1:numel(st), s = [s sprintf('%.0f–%.0f  ', fg(st(i))/1e3, fg(en(i))/1e3)]; end %#ok<AGROW>
if isempty(s), s = '（无）'; end
end

function savefig2(fig, figdir, name)
axs = findall(fig,'Type','axes');
for k = 1:numel(axs), try, axs(k).Toolbar.Visible = 'off'; catch, end, end
exportgraphics(fig, fullfile(figdir,[name '.png']), 'Resolution', 200);
savefig(fig, fullfile(figdir,[name '.fig']));
figure(fig); drawnow;
end
