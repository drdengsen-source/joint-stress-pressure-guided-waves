function C = T4_compare(S, T3, opts)
%T4_COMPARE  把 T4 的 FEM 频散（按支给出 f(k)）与 SAFE 参考解逐点比
%   C = T4_compare(S, T3)          S 来自 T4_sweep；T3 = load('T3_sweep.mat')
%   opts.fmin/fmax  只比这段频率（默认 20e3 / 800e3）
%   opts.tol        认为"配上了"的 c_p 相对差上限（默认 0.15）
%
% 做法：每一支 m 的 f_m(k) 单调，反过来插值得到 k(f_s)，于是同一频率下的
% c_p^FEM = 2πf_s/k，c_g^FEM、η_f^FEM 也在该支上插值。对每个 SAFE 点，
% 取 c_p 最接近的那一支作为对应支。
%
% 返回 C 是 table：f_Hz, cp_safe, cp_fem, dcp_pct, cg_safe, cg_fem, dcg_pct,
%                  eta_safe, eta_fem, branch, k_fem

if nargin < 3, opts = struct(); end
if ~isfield(opts,'fmin'), opts.fmin = 20e3; end
if ~isfield(opts,'fmax'), opts.fmax = 800e3; end
if ~isfield(opts,'tol'),  opts.tol  = 0.15; end

sel = T3.allF >= opts.fmin & T3.allF <= opts.fmax;
fs = T3.allF(sel);  cps = T3.allCP(sel);  cgs = T3.allCG(sel);  etas = T3.allETA(sel);
n = numel(fs);  nm = size(S.f,2);

CPf = nan(n,nm); CGf = nan(n,nm); ETf = nan(n,nm); Kf = nan(n,nm);
for m = 1:nm
    fm = S.f(:,m);  g = ~isnan(fm);
    if nnz(g) < 4, continue; end
    km = S.k(g);  fm = fm(g);
    [fm, iu] = unique(fm);  km = km(iu);
    if ~all(diff(fm) > 0), continue; end
    kk = interp1(fm, km, fs, 'pchip', NaN);
    Kf(:,m)  = kk;
    CPf(:,m) = 2*pi*fs./kk;
    CGf(:,m) = interp1(fm, S.cg(g,m), fs, 'pchip', NaN);
    ETf(:,m) = interp1(fm, S.eta(g,m), fs, 'pchip', NaN);
end
rel = abs(CPf./cps - 1);
[best, mb] = min(rel, [], 2, 'omitnan');
lin = sub2ind([n nm], (1:n).', mb);
ok = best <= opts.tol & ~isnan(best);

C = table(fs(ok), cps(ok), CPf(lin(ok)), 100*(CPf(lin(ok))./cps(ok)-1), ...
          cgs(ok), CGf(lin(ok)), 100*(CGf(lin(ok))./cgs(ok)-1), ...
          etas(ok), ETf(lin(ok)), mb(ok), Kf(lin(ok)), ...
    'VariableNames', {'f_Hz','cp_safe','cp_fem','dcp_pct','cg_safe','cg_fem','dcg_pct', ...
                      'eta_safe','eta_fem','branch','k_fem'});
C.Properties.UserData = struct('label', S.label, 'nmatched', nnz(ok), 'ntotal', n);
end
