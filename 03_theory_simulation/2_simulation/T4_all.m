function T4_all(port, kstep)
%T4_ALL  跑 T4 §1.2 的全部载荷组合，每个存 (k, f, c_p, c_g, eta_f) 与 §1.4 的 CSV
%   T4_all            默认端口 2038、k 步长 30 rad/m
%   T4_all(2038, 15)  更细的 k 网格
%
% 约定见 T4_build.m 头部。
% 输出：2_simulation\T4_runs\T4_all.mat（结构体 A，字段名 = 载荷编号）
%       2_simulation\T4_runs\T4_<id>.csv（列 f_Hz, k_ref_1perm, cp_nat, cg_nat, eta_f, load_id）
%
% 载荷表（比方案多两个：L2f / L2pf 是把闭端轴向力扣掉的版本，
% 因为 T3 的 S12 是"σz 作为独立未知量时对 p 的偏导"，必须在 σz = 0 下比）：
%   L0    基准                                L1    σz,ext = 100 MPa
%   L2    p = 20 MPa 闭端（σz = 50.4 自动）    L2p   同 L2 但液体参数冻结在 p=0（方案的 L2'）
%   L2f   p = 20 MPa 但 σz = 0                L2pf  同 L2f 且液体冻结
%   L3    100 MPa + 20 MPa（线性叠加检验）     L4    ΔT = +10 K
%   L5_0  水，基准                             L5    水，p = 20 MPa 闭端
%   L0w1  W1 线弹性基准                        L1w1  W1 + 100 MPa（分离三阶常数）
%   L0cf1400 / L0cf1500  改 c_f 看干净窗口位置随实测 c_f 怎么移

if nargin < 1 || isempty(port), port = 2038; end
if nargin < 2 || isempty(kstep), kstep = 30; end
here = fileparts(mfilename('fullpath'));  addpath(here); addpath(fullfile(here,'..','1_theory'));
outdir = fullfile(here,'T4_runs');  if ~exist(outdir,'dir'), mkdir(outdir); end
logdir = fullfile(here,'logs');     if ~exist(logdir,'dir'), mkdir(logdir); end
diary(fullfile(logdir,'T4_all.log'));
fprintf('==== T4_all port=%d kstep=%g  %s ====\n', port, kstep, char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')));
if isempty(which('mphstart')), addpath('D:\comsol64\Multiphysics\mli'); end
try, com.comsol.model.util.ModelUtil.tags; catch, mphstart(port); end

kg = 20:kstep:4220;
AB = 2.520833;                       % A_bore/A_wall
sz_open = sprintf('%.6g[MPa]', -20*AB);

% id, 模型键, sig_ext, p_in, frz, dT, c_f0覆盖
CS = { ...
 'L0',       'oilW2', '0[MPa]',  '0[MPa]',  '1', '0[K]',  ''; ...
 'L1',       'oilW2', '100[MPa]','0[MPa]',  '1', '0[K]',  ''; ...
 'L2',       'oilW2', '0[MPa]',  '20[MPa]', '1', '0[K]',  ''; ...
 'L2p',      'oilW2', '0[MPa]',  '20[MPa]', '0', '0[K]',  ''; ...
 'L2f',      'oilW2', sz_open,   '20[MPa]', '1', '0[K]',  ''; ...
 'L2pf',     'oilW2', sz_open,   '20[MPa]', '0', '0[K]',  ''; ...
 'L3',       'oilW2', '100[MPa]','20[MPa]', '1', '0[K]',  ''; ...
 'L4',       'oilW2', '0[MPa]',  '0[MPa]',  '1', '10[K]', ''; ...
 'L0cf1400', 'oilW2', '0[MPa]',  '0[MPa]',  '1', '0[K]',  '1400[m/s]'; ...
 'L0cf1500', 'oilW2', '0[MPa]',  '0[MPa]',  '1', '0[K]',  '1500[m/s]'; ...
 'L0w1',     'oilW1', '0[MPa]',  '0[MPa]',  '1', '0[K]',  ''; ...
 'L1w1',     'oilW1', '100[MPa]','0[MPa]',  '1', '0[K]',  ''; ...
 'L5_0',     'watW2', '0[MPa]',  '0[MPa]',  '1', '0[K]',  ''; ...
 'L5',       'watW2', '0[MPa]',  '20[MPa]', '1', '0[K]',  ''};

M = struct();
A = struct();
tAll = tic;
for i = 1:size(CS,1)
    id = CS{i,1};  key = CS{i,2};
    if ~isfield(M, key)
        switch key
            case 'oilW2', M.(key) = T4_build('wallmodel','W2','fluid','oil','quiet',true);
            case 'oilW1', M.(key) = T4_build('wallmodel','W1','fluid','oil','quiet',true);
            case 'watW2', M.(key) = T4_build('wallmodel','W2','fluid','water','quiet',true);
        end
        M.(key).study('std1').feature('eig').set('neigs', 14);
        fprintf('建模 %s\n', key);
    end
    mdl = M.(key);
    lp = struct('sig_ext',CS{i,3}, 'p_in',CS{i,4}, 'frz',CS{i,5}, 'dT',CS{i,6});
    if ~isempty(CS{i,7}), lp.c_f0 = CS{i,7}; else, mdl.param.set('c_f0', ternary(strcmp(key,'watW2'),'1482[m/s]','1450[m/s]')); end
    t0 = tic;
    S = T4_sweep(mdl, kg, lp, struct('label',id,'quiet',true));
    A.(id) = S;
    fprintf('%-9s 完成 %5.1f s  （%s，sig_ext=%s p=%s frz=%s dT=%s）\n', id, toc(t0), key, lp.sig_ext, lp.p_in, lp.frz, lp.dT);
    writecsv(fullfile(outdir, ['T4_' id '.csv']), S, id);
end
save(fullfile(outdir,'T4_all.mat'), 'A', 'kg', '-v7.3');
fprintf('==== 全部完成 %.1f min，存 %s ====\n', toc(tAll)/60, fullfile(outdir,'T4_all.mat'));
diary off
end

function writecsv(file, S, id)
% §1.4 固定列名：f_Hz, k_ref_1perm, cp_nat, cg_nat, eta_f, load_id
[nk, nm] = size(S.f);
K = repmat(S.k, 1, nm);
g = ~isnan(S.f) & S.f <= 900e3;
M = [S.f(g), K(g), S.cp(g), S.cg(g), S.eta(g)];
M = sortrows(M, [1 3]);
fid = fopen(file, 'w');
fprintf(fid, '%% T4 FEM (COMSOL 6.4 Floquet-Bloch 单胞)  载荷 %s  自然速度（k 按材料坐标）\n', id);
fprintf(fid, 'f_Hz,k_ref_1perm,cp_nat,cg_nat,eta_f,load_id\n');
for i = 1:size(M,1)
    fprintf(fid, '%.9g,%.9g,%.9g,%.9g,%.6g,%s\n', M(i,1), M(i,2), M(i,3), M(i,4), M(i,5), id);
end
fclose(fid);
end

function o = ternary(c,a,b), if c, o = a; else, o = b; end, end
