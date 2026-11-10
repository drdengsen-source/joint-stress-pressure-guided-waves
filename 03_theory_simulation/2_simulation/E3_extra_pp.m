function R = E3_extra_pp(datadir)
%E3_EXTRA_PP  瞬态仿真的三项补充分析，只读 E3_runs 里的 CSV，不连 COMSOL
%   R = E3_extra_pp                 默认 2_simulation\E3_runs
%
% Mo3  run A 的群速度用**频散补偿相位计时**核一遍：包络峰口径给出 c_g = 5208 m/s，比
%      SAFE 高 4.1%，原因是 5 周期突发在频散介质里包络峰前移和 200 ns 的输出步长。
%      这里在外壁测线上取两个相距 100 mm 的站（= 标距），由互谱相位直接量 k(f)，
%      再由 c_p = ω/k、c_g = dω/dk 给出与 SAFE 的对比，并给出用 SAFE 频散补偿后的残余延时。
% Mo4  ±z 对称性：run A / run B′ 各有一条 z_T 右侧测线（0.275–0.435 m）和一条镜像测线
%      （0.065–0.225 m）。比同一 |z − z_T| 上的包络峰幅值与到达时刻。
% M5   run B′（液组 280 kHz 5 周期窄带）：280 kHz 的 k 切片峰与能量份额，验证"280 kHz 处工作支占优"
%      这个结论不只依赖 run C 的宽带分解。
%
% 输出：<datadir>\E3_补充_结果.txt、<datadir>\E3_extra.mat（给出图脚本用），返回结构 R。
% 约定同 E3_postprocess.m：z 沿管轴，k > 0 沿 +z；z_T = 250 mm、z_R = 350 mm、标距 L = 100 mm。

if nargin < 1 || isempty(datadir)
    datadir = fullfile(fileparts(mfilename('fullpath')),'E3_runs');
end
here = fileparts(mfilename('fullpath'));
zT = 0.250;  L = 0.100;
T3 = load(fullfile(here,'..','1_theory','T3_sweep.mat'),'allF','allCP','allCG','allETA');
T3.k = 2*pi*T3.allF./T3.allCP;

fid = fopen(fullfile(datadir,'E3_补充_结果.txt'),'w','n','UTF-8');
pr = @(varargin) fprintf(fid, varargin{:});
pr('E3 补充分析  %s\n', char(datetime('now','Format','yyyy-MM-dd HH:mm')));
pr('数据目录 %s；z_T = %.0f mm，标距 L = %.0f mm。这是仿真结果，不是实验数据。\n\n', datadir, zT*1e3, L*1e3);
R = struct();

%% ================= Mo3：run A 的相位计时 c_g =================
pr('== Mo3  run A（壁组 520 kHz）：频散补偿相位计时核群速度 ==\n');
[zA, tA, UA] = readLine(fullfile(datadir,'runA_line_uz.csv'));
f0 = 520e3;
% SAFE 壁传族在 460–580 kHz 的 k(f)：用于补偿与对照
wsel = T3.allETA < 0.5 & T3.allF >= 3e5 & T3.allF <= 8e5;
fS = T3.allF(wsel);  kS = T3.k(wsel);  [fS,o] = sort(fS);  kS = kS(o);
% 两个站：z1、z2 相距 L；避开发射近场（z1 = z_T + 35 mm）
z1 = zT + 0.035;  z2 = z1 + L;
[~,i1] = min(abs(zA - z1));  [~,i2] = min(abs(zA - z2));
dz = zA(i2) - zA(i1);
s1 = UA(i1,:).';  s2 = UA(i2,:).';
dt = tA(2) - tA(1);  nt = numel(tA);  nfft = 2^nextpow2(4*nt);
fax = (0:nfft-1).'/(nfft*dt);
S1 = fft(s1, nfft);  S2 = fft(s2, nfft);
band = fax >= 470e3 & fax <= 570e3;
fb = fax(band);
% 约定：MATLAB 的 fft 核是 e^{-i2πft}，所以 cos(kz − ωt) 在**正频率**上的系数 ∝ e^{−ikz}
%（与 E3_postprocess.m 的 fk_map 里 k = −2π·κ 同一个约定）→ 互谱相位 = −k·Δz。
X  = S2(band).*conj(S1(band));
ph = unwrap(angle(X));
kSAFE = interp1(fS, kS, fb, 'pchip');
n2pi  = round(mean(kSAFE.*dz + ph)/(2*pi));         % unwrap 只给相对相位，绝对支由 SAFE 定
kmeas = (2*pi*n2pi - ph)/dz;
wgt   = abs(S2(band)).*abs(S1(band));  wgt = wgt/max(wgt);
% c_p、c_g：k(ω) 的加权二次拟合（自变量中心化并归一，否则 ω² ~ 10^13 会让最小二乘秩亏）
om0 = 2*pi*f0;
fitk = @(sel) wfit((2*pi*fb(sel) - om0)/om0, kmeas(sel), wgt(sel));
p100 = fitk(true(size(fb)));
cp_meas = om0/p100(1);
cp_safe = interp1(fS, 2*pi*fS./kS, f0, 'pchip');
cg_safe = interp1(T3.allF(wsel), T3.allCG(wsel), f0, 'pchip');
% c_g 随拟合带宽（支的曲率 + 突发带宽有限 → 宽带拟合是带内平均，不是 520 kHz 的局部导数）
bws = [20 30 40 60 100]*1e3;  cg_bw = nan(size(bws));
for i = 1:numel(bws)
    sel = abs(fb - f0) <= bws(i)/2;
    pp = fitk(sel);  cg_bw(i) = om0/pp(2);
end
cg_meas = cg_bw(bws == 40e3);
% 频散补偿互相关：把 z1 的解析信号用 SAFE 的 k(f) 推到 z2（正频率乘 e^{−ik·Δz}），与 z2 实测比
kf = zeros(nfft,1);  pos = fax < 1/(2*dt);
gp = pos & fax >= min(fS) & fax <= max(fS);
kf(gp) = interp1(fS, kS, fax(gp), 'pchip');
A1 = S1;  A1(~pos) = 0;  A1(pos) = 2*A1(pos);  A1(1) = A1(1)/2;    % 解析信号谱
A2 = S2;  A2(~pos) = 0;  A2(pos) = 2*A2(pos);  A2(1) = A2(1)/2;
a1p = ifft(A1.*exp(-1i*kf*dz));  a1p = a1p(1:nt);
a2a = ifft(A2);  a2a = a2a(1:nt);
[lag_ns, rho, ph_res] = xcorrpeak(a1p, a2a, dt);
% 不补偿的包络峰口径（+4.1% 的那个量）
e1 = envelope_fft(s1);  e2 = envelope_fft(s2);
[~,j1] = max(e1);  [~,j2] = max(e2);
cg_env = dz/((j2-j1)*dt);
pr('  两站 z1 = %.1f mm、z2 = %.1f mm（Δz = %.1f mm = 标距）；互谱相位给 k(f)，权重 = |S_1||S_2|\n', zA(i1)*1e3, zA(i2)*1e3, dz*1e3);
pr('    f/kHz    k 实测      k SAFE     c_p 实测   c_p SAFE    差      权重\n');
for f = [480 500 520 540 560]*1e3
    [~,j] = min(abs(fb - f));
    pr('    %5.0f  %.4f rad/mm %.4f    %7.1f  %7.1f  %+6.2f%%   %.2f\n', f/1e3, kmeas(j)/1e3, kSAFE(j)/1e3, ...
        2*pi*fb(j)/kmeas(j), 2*pi*fb(j)/kSAFE(j), 100*(kSAFE(j)/kmeas(j)-1), wgt(j));
end
pr('  520 kHz：c_p = %.1f m/s（SAFE %.1f，%+.2f%%）\n', cp_meas, cp_safe, 100*(cp_meas/cp_safe-1));
pr('  c_g = dω/dk，随拟合带宽：%s（SAFE 的能流口径 c_g = P/E 给 %.0f m/s）\n', ...
    strjoin(arrayfun(@(i) sprintf('%.0f kHz → %.0f', bws(i)/1e3, cg_bw(i)), 1:numel(bws), 'UniformOutput', false), '，'), cg_safe);
pr('    （SAFE 自己的 k(f) 样点在 ±30 kHz 上差分只给 4835 m/s，所以宽带拟合的偏低是支的曲率，不是模型差）\n');
pr('  同一对站不补偿的包络峰口径：c_g = %.0f m/s（%+.1f%%）\n', cg_env, 100*(cg_env/cg_safe-1));
pr('  **频散补偿互相关**（用 SAFE 的 k(f) 把 z1 推到 z2）：残余群延时 %+.0f ns，占模型群延时 %.2f µs 的 %+.2f%%；\n', ...
    lag_ns, dz/cg_safe*1e6, 100*lag_ns*1e-9/(dz/cg_safe));
pr('    相关系数 %.4f，峰处相位残差 %.3f rad = %+.0f ns（相延时 %.2f µs 的 %+.2f%%）。\n', ...
    rho, ph_res, ph_res/(2*pi*f0)*1e9, dz/cp_safe*1e6, 100*(ph_res/(2*pi*f0))/(dz/cp_safe));
pr('    结论：包络峰口径的 +4~6%% 是 5 周期突发在频散支上的包络前移，不是模型偏差；\n');
pr('    §2.4 规定的频散补偿相位计时把群延时还原到 %.1f%% 以内。\n\n', abs(100*lag_ns*1e-9/(dz/cg_safe)));
R.Mo3 = struct('z1',zA(i1),'z2',zA(i2),'dz',dz,'f',fb,'k',kmeas,'kSAFE',kSAFE,'w',wgt, ...
    'cp_meas',cp_meas,'cp_safe',cp_safe,'cg_bw',cg_bw,'bws',bws,'cg_meas',cg_meas,'cg_safe',cg_safe, ...
    'cg_env',cg_env,'lag_ns',lag_ns,'rho',rho,'ph_res',ph_res);

%% ================= Mo4：±z 对称性 =================
pr('== Mo4  ±z 对称性（发射力分布关于 z_T 对称 → 两个方向幅值应相同）==\n');
R.Mo4 = struct();
for cs = {{'A','uz','壁组 520 kHz，外壁 u_z'}, {'BP','ur','液组 280 kHz，外壁 u_r'}}
    id = cs{1}{1};  comp = cs{1}{2};  nm = cs{1}{3};
    fp = fullfile(datadir, sprintf('run%s_line_%s.csv', id, comp));
    fn = fullfile(datadir, sprintf('run%s_lineneg_%s.csv', id, comp));
    if ~exist(fp,'file') || ~exist(fn,'file')
        pr('  run %s：缺 %s 或镜像测线，跳过\n', id, comp);  continue;
    end
    [zp, tp, Up] = readLine(fp);
    [zn, ~,  Un] = readLine(fn);
    dpos = zp - zT;  dneg = zT - zn;                       % 到发射中心的距离
    dsel = 0.025:0.005:0.085;                              % 只比 25–85 mm：避开发射近场和接收探头（z_R−z_T = 100 mm）
    ap = nan(size(dsel));  an = ap;  tpk_p = ap;  tpk_n = ap;
    for i = 1:numel(dsel)
        [~,jp] = min(abs(dpos - dsel(i)));  [~,jn] = min(abs(dneg - dsel(i)));
        ep = envelope_fft(Up(jp,:).');  en = envelope_fft(Un(jn,:).');
        [ap(i),kp] = max(ep);  [an(i),kn] = max(en);
        tpk_p(i) = tp(kp);  tpk_n(i) = tp(kn);
    end
    rat = an./ap;
    pr('  run %s（%s）：|z−z_T| = %.0f–%.0f mm 上 %d 个站\n', id, nm, dsel(1)*1e3, dsel(end)*1e3, numel(dsel));
    pr('    包络峰幅值比 −z/+z：均值 %.4f，范围 %.4f–%.4f（%.2f dB）\n', mean(rat), min(rat), max(rat), 20*log10(mean(rat)));
    pr('    包络峰到达时刻差 −z 减 +z：均值 %+.0f ns，最大 %+.0f ns（输出步长 %.0f ns）\n', ...
        mean(tpk_n-tpk_p)*1e9, max(abs(tpk_n-tpk_p))*1e9, (tp(2)-tp(1))*1e9);
    R.Mo4.(['run' id]) = struct('d',dsel,'a_pos',ap,'a_neg',an,'ratio',rat,'dt',tpk_n-tpk_p);
end
pr('\n');

%% ================= M5：run B′（280 kHz 窄带）=================
fgBP = fullfile(datadir,'runBP_global.csv');
if exist(fgBP,'file')
    pr('== M5  run B′（液组 280 kHz，5 周期）：280 kHz 的 k 切片份额 ==\n');
    G = readmatrix(fgBP,'FileType','text','CommentStyle','%');
    t = G(:,1);  Itx = G(:,2);  Vw = G(:,3);  Vl = G(:,4);
    envI = envelope_fft(Itx);  [~,i0] = max(envI);  t0 = t(i0);
    envL = envelope_fft(Vl);  envW = envelope_fft(Vw);
    pr('  记录 %.0f µs，步长 %.0f ns；发射突发包络峰 t0 = %.1f µs；发射电流峰值 %.2f A\n', ...
        t(end)*1e6, (t(2)-t(1))*1e9, t0*1e6, max(abs(Itx)));
    pr('  串扰：壁组通道最大包络 / 液组通道最大包络 = %.2f（%.1f dB）\n', max(envW)/max(envL), 20*log10(max(envW)/max(envL)));
    % 时间门 60–86 µs 内的包络峰
    gate = t >= t0 + 55e-6 & t <= t0 + 80e-6;
    [~,jg] = max(envL.*gate);
    pr('  时间门（t0 + 55–80 µs）内包络峰在 %.1f µs（到达 %.1f µs → 表观 c_g = %.0f m/s；\n', ...
        t(jg)*1e6, (t(jg)-t0)*1e6, L/(t(jg)-t0));
    pr('    SAFE 280 kHz：工作支 c_g 1483 → 到达 %.1f µs；上支 1197 → %.1f µs）\n', L/1483*1e6, L/1197*1e6);
    % k 切片
    [zl, tl, U] = readLine(fullfile(datadir,'runBP_line_ur.csv'));
    [FK, fk_f, fk_k] = fk_map(U, zl, tl);
    f0 = 280e3;
    fsel = abs(fk_f - f0) <= 10e3;
    sl = mean(FK(fsel,:),1);
    kpos = fk_k > 2*pi*f0/2500 & fk_k < 2*pi*f0/1100;      % c_p 在 1100–2500 m/s
    [pkk,ikk] = findpeaks_simple(sl(:).*kpos(:), 3, 0.1*max(sl(:).*kpos(:)));
    cps = 2*pi*f0./fk_k(ikk);
    share = pkk.^2/sum(pkk.^2);
    cand = [1979.1 1565.7 1339.8];  cgc = [1197.4 1482.9 1722.8];  cnm = {'上支','工作支','下支'};
    pr('  280 kHz ± 10 kHz 的 k 切片峰：\n');
    for j = 1:numel(cps)
        [~,im] = min(abs(cand - cps(j)));
        pr('    c_p = %6.0f m/s（k = %.3f rad/mm），能量份额 %.2f → SAFE 的%s（%.0f m/s，差 %+.1f%%，c_g %.0f → 到达 %.1f µs）\n', ...
            cps(j), fk_k(ikk(j))/1e3, share(j), cnm{im}, cand(im), 100*(cps(j)/cand(im)-1), cgc(im), L/cgc(im)*1e6);
    end
    % 沿三支脊线积分的份额（与 run C 同口径，但只在 275–285 kHz）
    br = branches280();
    fs2 = fk_f >= 275e3 & fk_f <= 285e3;  fB = fk_f(fs2);  A2 = FK(fs2,:).^2;
    Em = zeros(numel(fB), numel(br));
    for m = 1:numel(br)
        cpB = interp1(br(m).f, br(m).cp, fB, 'pchip', 'extrap');
        kB = 2*pi*fB./cpB;
        for i = 1:numel(fB)
            Em(i,m) = sum(A2(i, abs(fk_k - kB(i)) <= 0.10*kB(i)));
        end
    end
    shb = sum(Em,1)/sum(Em,'all');
    pr('  沿 SAFE 三支脊线 ±10%% 积分（275–285 kHz，与 run C 同口径）：%s\n', ...
        strjoin(arrayfun(@(m) sprintf('%s %.2f', br(m).name, shb(m)), 1:numel(br), 'UniformOutput', false), '，'));
    pr('  对照 run B（300 kHz 窄带）：上支 0.44 / 工作支 0.43 / 下支 0.11；run C（宽带 250–400 kHz 积分）：0.68 / 0.16 / 0.15\n\n');
    R.M5 = struct('t',t,'Itx',Itx,'Vl',Vl,'Vw',Vw,'t0',t0,'env',envL, ...
        'k',fk_k,'slice',sl,'cps',cps,'share',share,'shb',shb,'brnames',{{br.name}}, ...
        'fk',struct('f',fk_f,'k',fk_k,'A',FK));
else
    pr('== M5  没找到 %s，跳过 ==\n\n', fgBP);
end

fclose(fid);
save(fullfile(datadir,'E3_extra.mat'),'R','-v7.3');
fprintf('E3 补充分析完成：%s\n', fullfile(datadir,'E3_补充_结果.txt'));
end

%% ---------- 局部函数 ----------
function br = branches280()
br(1) = struct('name','上支',   'f',[250 280 300 350]*1e3, 'cp',[2181.4 1979.1 1904.1 1810.8]);
br(2) = struct('name','工作支', 'f',[250 280 300 350]*1e3, 'cp',[1578.6 1565.7 1560.4 1550.5]);
br(3) = struct('name','下支',   'f',[250 280 300 350]*1e3, 'cp',[1304.3 1339.8 1359.4 1395.1]);
end

function w = tukey(n, r)
w = ones(n,1);  m = max(2, round(r*n));
w(1:m) = 0.5*(1-cos(pi*(0:m-1).'/(m-1)));  w(end-m+1:end) = flipud(w(1:m));
end

function p = wfit(x, y, w)
% k(ω) ≈ p(1) + p(2)·x + p(3)·x²，x = (ω−ω0)/ω0；返回时把 p(2) 换算成 dk/dω 的倒数用
x = x(:);  y = y(:);  W = sqrt(w(:)/max(w));
V = [ones(size(x)) x x.^2];
p = (W.*V)\(W.*y);
end

function [lag_ns, rho, ph_res] = xcorrpeak(a, b, dt)
% a、b 是解析信号：|互相关| 的峰给群延时残差，峰处相位给相延时残差
a = a(:);  b = b(:);  n = numel(a);  nf = 2^nextpow2(2*n);
c = ifft(fft(a,nf).*conj(fft(b,nf)));
c = [c(nf-n+2:nf); c(1:n)];                 % lag = −(n−1) … (n−1)
lags = (-(n-1):(n-1)).';
y = abs(c);  [~, im] = max(y);
if im > 1 && im < numel(y)
    d = (y(im-1) - y(im+1))/(2*(y(im-1) - 2*y(im) + y(im+1)));
else
    d = 0;
end
lag_ns = (lags(im) + d)*dt*1e9;
rho = y(im)/sqrt(sum(abs(a).^2)*sum(abs(b).^2));
ph_res = angle(c(im));
end

function [z, t, U] = readLine(file)
txt = fileread(file);
lines = regexp(txt, '\r?\n', 'split');
hdr = lines(startsWith(strtrim(lines), '%'));
t = [];
for i = 1:numel(hdr)
    tok = regexp(hdr{i}, 't=([-+0-9.eE]+)', 'tokens');
    if ~isempty(tok), t = cellfun(@(c) str2double(c{1}), tok); break; end
end
M = readmatrix(file, 'FileType','text', 'CommentStyle','%');
z = M(:,2);  U = M(:,3:end);
if numel(t) ~= size(U,2), error('readLine: 时刻数 %d 与数据列数 %d 不符', numel(t), size(U,2)); end
[z, iz] = sort(z);  U = U(iz,:);  t = t(:);
end

function env = envelope_fft(x)
x = x(:);  n = numel(x);  X = fft(x);  h = zeros(n,1);
if mod(n,2)==0, h([1 n/2+1])=1; h(2:n/2)=2; else, h(1)=1; h(2:(n+1)/2)=2; end
env = abs(ifft(X.*h));
end

function [pk, ipk] = findpeaks_simple(x, mindist, minh)
x = x(:);  cand = find(x(2:end-1) > x(1:end-2) & x(2:end-1) >= x(3:end) & x(2:end-1) >= minh) + 1;
[~, ord] = sort(x(cand), 'descend');  cand = cand(ord);  keep = [];
for i = 1:numel(cand)
    if all(abs(cand(i) - keep) >= mindist), keep(end+1) = cand(i); end %#ok<AGROW>
end
ipk = sort(keep);  pk = x(ipk);  ipk = ipk(:);  pk = pk(:);
if isempty(ipk), [pk, ipk] = max(x); end
end

function [A, f, k] = fk_map(U, z, t)
if size(U,2) > 3000, U = U(:,1:2:end); t = t(1:2:end); end
[nz, nt] = size(U);
wz = 0.5*(1 - cos(2*pi*(0:nz-1).'/(nz-1)));
wt = ones(1,nt);  nr = max(2, round(0.05*nt));
wt(1:nr) = 0.5*(1 - cos(pi*(0:nr-1)/(nr-1)));  wt(end-nr+1:end) = fliplr(wt(1:nr));
Uw = (U - mean(U,2)).*wz.*wt;
nk = 2^nextpow2(4*nz);  nf = 2^nextpow2(2*nt);
F = fft2(Uw, nk, nf);
dz = z(2)-z(1);  dt = t(2)-t(1);
kax = (0:nk-1)/(nk*dz);  kax(kax >= 1/(2*dz)) = kax(kax >= 1/(2*dz)) - 1/dz;
fax = (0:nf-1)/(nf*dt);
fsel = fax <= 1/(2*dt);  f = fax(fsel).';
A = abs(F(:,fsel)).';
k = -2*pi*kax;
[k, ik] = sort(k);  A = A(:,ik);
end
