function S = T4_sweep(model, kgrid, loadpar, opts)
%T4_SWEEP  对 T4 单胞做 Floquet 波数扫描，返回按支（列指标）组织的频散
%   S = T4_sweep(model, kgrid)
%   S = T4_sweep(model, kgrid, struct('p_in','20[MPa]'), struct('label','L2'))
%
%   S.k (nk×1)；S.f / S.cp / S.cg / S.eta (nk×nm)，第 m 列是第 m 支。
%   固定 k 时把本征值按频率升序排，得到的列指标就是一支连续频散曲线：
%   耦合系统只有 avoided crossing，本征值不会真交叉，所以不需要额外的模态跟踪。
%   c_g = dω/dk 直接在 k 网格上做中心差分（比再解一次 k+δk 快一倍，网格够密时更稳）。
%
%   opts.label / opts.quiet
%
% 约定（方案 §0.1）：Floquet 相位按材料坐标 → k 是参考波数 → c_p = 2πf/k 是**自然速度**。
% eta_f 与 SAFE 同口径 (Vf+Tf)/E；无耗本征模均分 Ts+Tf = Vs+Vf ⇒ E = 2(Ts+Tf)。

if nargin < 3 || isempty(loadpar), loadpar = struct(); end
if nargin < 4, opts = struct(); end
if ~isfield(opts,'label'), opts.label = 'L?'; end
if ~isfield(opts,'quiet'), opts.quiet = false; end

fn = fieldnames(loadpar);
for i = 1:numel(fn), model.param.set(fn{i}, loadpar.(fn{i})); end

eTs = '0.25*rho_alT*(2*pi*freq)^2*(abs(u)^2+abs(w)^2)';        % 管壁动能密度（时均）
eTf = '0.25*(abs(pr)^2+abs(pz)^2)/(rho_fT*(2*pi*freq)^2)';     % 液体动能密度
eVf = '0.25*abs(p)^2/(rho_fT*c_fT^2)';                         % 液体势能密度

kgrid = sort(kgrid(:));  nk = numel(kgrid);
F = []; ETA = [];
for i = 1:nk
    model.param.set('kz', sprintf('%.12g[1/m]', kgrid(i)));
    if i == 1, model.sol('sol1').runAll;         % 第一次要算静力步（预应力）
    else,      model.sol('sol1').runFrom('st2'); % 之后只重解本征值：kz 不影响静力步
    end
    f = real(mphglobal(model,'freq','dataset','dset1','solnum','all'));  f = f(:);
    Ts = rcol(mphint2(model, eTs, 2, 'dataset','dset1','selection','sel_wall', 'intvolume','on','solnum','all'));
    [tf, vf] = mphint2(model, {eTf, eVf}, 2, 'dataset','dset1','selection','sel_fluid','intvolume','on','solnum','all');
    Tf = rcol(tf);  Vf = rcol(vf);
    e = (Vf + Tf)./(2*(Ts + Tf));
    [f, o] = sort(f);  e = e(o);
    if isempty(F), nm = numel(f);  F = nan(nk,nm);  ETA = nan(nk,nm); end
    n = min(numel(f), nm);
    F(i,1:n) = f(1:n);  ETA(i,1:n) = e(1:n);
    if ~opts.quiet && (mod(i,25)==0 || i==nk)
        fprintf('  %s: %d/%d  k = %.0f rad/m，最低支 %.1f kHz\n', opts.label, i, nk, kgrid(i), F(i,1)/1e3);
    end
end
CG = nan(size(F));
for m = 1:size(F,2)
    g = ~isnan(F(:,m));
    if nnz(g) > 2, CG(g,m) = 2*pi*gradient(F(g,m), kgrid(g)); end
end
S = struct('label',opts.label,'k',kgrid,'f',F,'cp',2*pi*F./kgrid,'cg',CG,'eta',ETA,'loadpar',loadpar);
end

function y = rcol(x), y = real(x(:)); end
