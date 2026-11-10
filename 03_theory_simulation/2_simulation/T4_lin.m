function T4_lin(port, kstep)
%T4_LIN  载荷线性度：在 p = 5 / 10 / 20 MPa（σz = 0）和 σz = 50 / 100 MPa（p = 0）各做一次 k 扫描，
%        用固定频率口径给出 K 随载荷幅度的变化 —— 判"标定曲线在量程内线性吗"。
%        方案 §1.3 的 L3 只查交叉项叠加，这一项查单载荷本身的线性度。
% 输出：2_simulation\T4_runs\T4_lin.mat + 追加到 T4_结果.txt
if nargin < 1 || isempty(port), port = 2038; end
if nargin < 2 || isempty(kstep), kstep = 30; end
here = fileparts(mfilename('fullpath')); addpath(here); addpath(fullfile(here,'..','1_theory'));
outdir = fullfile(here,'T4_runs');
diary(fullfile(here,'logs','T4_lin.log'));
if isempty(which('mphstart')), addpath('D:\comsol64\Multiphysics\mli'); end
try, com.comsol.model.util.ModelUtil.tags; catch, mphstart(port); end
kg = 20:kstep:4220;  AB = 2.520833;
M = T4_build('quiet',true);  M.study('std1').feature('eig').set('neigs',14);
B = struct();
B.base = T4_sweep(M, kg, struct('sig_ext','0[MPa]','p_in','0[MPa]','frz','1','dT','0[K]'), struct('label','base','quiet',true));
for p = [5 10]
    id = sprintf('p%d', p);
    B.(id) = T4_sweep(M, kg, struct('sig_ext',sprintf('%.8g[MPa]',-p*AB),'p_in',sprintf('%d[MPa]',p),'frz','1','dT','0[K]'), ...
                      struct('label',id,'quiet',true));
    fprintf('%s 完成\n', id);
end
B.s50 = T4_sweep(M, kg, struct('sig_ext','50[MPa]','p_in','0[MPa]','frz','1','dT','0[K]'), struct('label','s50','quiet',true));
fprintf('s50 完成\n');
save(fullfile(outdir,'T4_lin.mat'),'B','kg','-v7.3');
fprintf('存 %s\n', fullfile(outdir,'T4_lin.mat'));
diary off
end
