function model = T4_build(varargin)
%T4_BUILD  建 T4 的 Floquet–Bloch 单胞模型（2D 轴对称：solid 管壁 + acpr 液体）
%   model = T4_build()                     默认 W2（Murnaghan 超弹）、油、壁 8 层 / 液 40 层
%   model = T4_build('wallmodel','W1')     线弹性 + 几何非线性（分离三阶常数用）
%   model = T4_build('nw',16)              网格收敛用
%   model = T4_build('fluid','water')      L5
%
% 参数：1_theory\params_pipe.m。
% 单位 SI。坐标 r 径向、z 轴向。应力拉为正。
%
% 结构（对应方案 §0）：
%   Study 1 步 1  Stationary（几何非线性开）：z=0 对称、z=L 轴向牵引 σz,ext + p·A_bore/A_wall、
%                 内壁随动压力 p。Floquet 周期条件在这一步关掉（它与"单胞被拉长"不相容）。
%   Study 1 步 2  Eigenfrequency（以步 1 为线性化点）：打开 solid/acpr 的 Floquet 周期条件
%                 （相位 kz·L 按材料坐标 → 给出**自然速度** c_nat，见方案 §0.1）；
%                 关掉 z=0 对称与 z=L 牵引（它们与周期条件冲突/不属于波问题）；
%                 **保留内壁随动压力**——随动压力刚度正是方案要量的那一项。
%   液体域固定在参考孔径（不做动网格），c_f(p,T)、ρ_f(p,T) 用物理值（方案 §0.3）。
%
% 载荷由参数控制（T4_sweep 负责设）：sig_ext、p_in、dT、frz（=0 冻结液体参数 = L2'）。

opt = struct('wallmodel','W2', 'fluid','oil', 'nw',8, 'nf',40, 'nz',4, 'Lcell',0.5e-3, 'quiet',false, ...
             'ri',[], 'ro',[]);
for k = 1:2:numel(varargin), opt.(varargin{k}) = varargin{k+1}; end
say = @(varargin) fprintf([varargin{1} '\n'], varargin{2:end});
if opt.quiet, say = @(varargin) []; end

Pp = params_pipe('7075', opt.fluid);
G = Pp.geom;  W = Pp.wall;  F = Pp.fluid;
% 几何更新检查（方案 §1.3 最后两行）：直接改参考孔径/外径，不用动网格
if ~isempty(opt.ri), G.r_i = opt.ri; end
if ~isempty(opt.ro), G.r_o = opt.ro; end

import com.comsol.model.util.*
model = ModelUtil.create(['T4_' char(java.util.UUID.randomUUID.toString.substring(0,6))]);
model.modelNode.create('comp1');

%% ---------- 参数 ----------
p = { ...
 'r_i',   sprintf('%.10g[m]',G.r_i),   '内壁半径'; ...
 'r_o',   sprintf('%.10g[m]',G.r_o),   '外壁半径'; ...
 'L_c',   sprintf('%.10g[m]',opt.Lcell),'单胞长（Floquet 相位 kz*L_c）'; ...
 'AB_AW', sprintf('%.10g',G.A_bore/G.A_wall), '闭端换算 A_bore/A_wall'; ...
 'lam0',  sprintf('%.10g[Pa]',W.lam),  'Lame lambda（Stobbe 2005）'; ...
 'mu0',   sprintf('%.10g[Pa]',W.mu),   'Lame mu'; ...
 'lM',    sprintf('%.10g[Pa]',W.l),    'Murnaghan l'; ...
 'mM',    sprintf('%.10g[Pa]',W.m),    'Murnaghan m'; ...
 'nM',    sprintf('%.10g[Pa]',W.n),    'Murnaghan n'; ...
 'rho_al0',sprintf('%.10g[kg/m^3]',W.rho), '7075 密度'; ...
 'alph_al',sprintf('%.10g[1/K]',W.alpha),  '线膨胀系数'; ...
 'dlnE_dT',sprintf('%.10g[1/K]',W.dlnE_dT),'(1/E)dE/dT'; ...
 'c_f0',  sprintf('%.10g[m/s]',F.c),   '液体声速 @p0,T0'; ...
 'rho_f0',sprintf('%.10g[kg/m^3]',F.rho), '液体密度 @p0,T0'; ...
 'a_p',   sprintf('%.10g[1/Pa]',F.a_p),'(1/c)dc/dp 等温（Beyer 修正）'; ...
 'K_T',   sprintf('%.10g[Pa]',F.KT),   '液体等温体积模量'; ...
 'dcf_dT',sprintf('%.10g[m/(s*K)]',F.dc_dT),'dc_f/dT'; ...
 'beta_f',sprintf('%.10g[1/K]',F.beta),'液体体膨胀系数'; ...
 'sig_ext','0[MPa]',  '外加轴向应力 F/A_wall'; ...
 'p_in',  '0[MPa]',   '内压'; ...
 'dT',    '0[K]',     '温升'; ...
 'frz',   '1',        '=0 把液体参数冻结在 p=0（L2''）'; ...
 'kz',    '500[1/m]', 'Floquet 轴向波数（按参考长度）'; ...
 'sigz_end','sig_ext + p_in*AB_AW', '端面轴向牵引（含闭端项）'};
for i = 1:size(p,1), model.param.set(p{i,1}, p{i,2}, p{i,3}); end

% 随温度/压力变化的材料量（全局变量）
vg = model.variable.create('vg');  vg.model('comp1');
vg.set('E_alT',  'mu0*(3*lam0+2*mu0)/(lam0+mu0)*(1+dlnE_dT*dT)', '杨氏模量(T)');
vg.set('nu_al',  'lam0/(2*(lam0+mu0))',                          '泊松比');
vg.set('rho_alT','rho_al0*(1-3*alph_al*dT)',                     '密度(T)');
vg.set('lamT',   'lam0*(1+dlnE_dT*dT)',                          'lambda(T)');
vg.set('muT',    'mu0*(1+dlnE_dT*dT)',                           'mu(T)');
vg.set('c_fT',   '(c_f0+dcf_dT*dT)*(1+a_p*p_in*frz)',            '液体声速(p,T)');
vg.set('rho_fT', 'rho_f0*(1-beta_f*dT)*(1+p_in*frz/K_T)',        '液体密度(p,T)');

%% ---------- 几何 ----------
g = model.geom.create('geom1',2);  g.axisymmetric(true);  g.lengthUnit('m');
g.create('rf','Rectangle');  g.feature('rf').set('pos',[0 0]);      g.feature('rf').set('size',[G.r_i opt.Lcell]);
g.create('rw','Rectangle');  g.feature('rw').set('pos',[G.r_i 0]);  g.feature('rw').set('size',[G.r_o-G.r_i opt.Lcell]);
g.run;
e = 1e-9;  L = opt.Lcell;
dFl = mphselectbox(model,'geom1',[-e G.r_i+e; -e L+e],'domain');
dWa = mphselectbox(model,'geom1',[G.r_i-e G.r_o+e; -e L+e],'domain');
bFl0 = mphselectbox(model,'geom1',[-e G.r_i+e; -e e],'boundary');       % 液体 z=0
bFlL = mphselectbox(model,'geom1',[-e G.r_i+e; L-e L+e],'boundary');    % 液体 z=L
bWa0 = mphselectbox(model,'geom1',[G.r_i-e G.r_o+e; -e e],'boundary');
bWaL = mphselectbox(model,'geom1',[G.r_i-e G.r_o+e; L-e L+e],'boundary');
bFSI = mphselectbox(model,'geom1',[G.r_i-e G.r_i+e; -e L+e],'boundary');
say('几何：液体域 %s，管壁域 %s；界面边界 %s', mat2str(dFl), mat2str(dWa), mat2str(bFSI));
sl = model.selection.create('sel_fluid','Explicit'); sl.geom('geom1',2); sl.set(dFl); sl.label('液体域');
sl = model.selection.create('sel_wall','Explicit');  sl.geom('geom1',2); sl.set(dWa); sl.label('管壁域');

%% ---------- 材料 ----------
iso = @(x) {x,'0','0','0',x,'0','0','0',x};
mw = model.material.create('mat_w','Common');  mw.selection.set(dWa);  mw.label('7075-T6');
mw.propertyGroup('def').set('density','rho_alT');
mw.propertyGroup.create('Enu','Young''s modulus and Poisson''s ratio');
mw.propertyGroup('Enu').set('E','E_alT');  mw.propertyGroup('Enu').set('nu','nu_al');
mf = model.material.create('mat_f','Common');  mf.selection.set(dFl);  mf.label(F.name);
mf.propertyGroup('def').set('density','rho_fT');
mf.propertyGroup('def').set('soundspeed','c_fT');

%% ---------- solid ----------
sd = model.physics.create('solid','SolidMechanics','geom1');
sd.selection.set(dWa);
% 材料模型：W1 线弹性（默认 lemm1）或 W2 Murnaghan 超弹
if strcmpi(opt.wallmodel,'W2')
    % 默认线弹性节点不能删；把超弹节点放在它后面即可覆盖同一域
    hm = sd.create('hmm1','HyperelasticModel',2);
    hm.selection.set(dWa);
    hm.set('MaterialModel','Murnaghan');
    try, hm.set('IsotropicOption','Lame'); catch, end
    hm.set('lambLame_mat','userdef');  hm.set('lambLame','lamT');
    hm.set('muLame_mat','userdef');    hm.set('muLame','muT');
    hm.set('l_mat','userdef');         hm.set('l','lM');
    hm.set('m_mat','userdef');         hm.set('m','mM');
    hm.set('n_mat','userdef');         hm.set('n','nM');
    hm.label('W2 Murnaghan 超弹');
else
    sd.feature('lemm1').label('W1 线弹性 + 几何非线性');
end
% 静力步：z=0 对称、z=L 轴向牵引
sym0 = sd.create('symz0','SymmetrySolid',1);  sym0.selection.set(bWa0);  sym0.label('静力步 z=0 对称');
blz  = sd.create('blz','BoundaryLoad',1);     blz.selection.set(bWaL);
blz.set('forceType','ForceArea');             % 按参考面积（σz,ext = F/A_wall 就是这个口径）
% 注意：forceType = ForceArea 时真正生效的向量属性是 forceReferenceArea，不是 F（F 是残留槽位，
% 写进去读得回来但不参与计算 —— 踩过一次：静力步位移恒为 0）
blz.set('forceReferenceArea_src','userdef');
blz.set('forceReferenceArea',{'0';'0';'sigz_end'});
blz.label('静力步 端面轴向牵引');
% 内壁随动压力（两步都开 —— 随动压力刚度是方案要量的一项）
blp = sd.create('blp','BoundaryLoad',1);      blp.selection.set(bFSI);
blp.set('forceType','FollowerPressure');      blp.set('pressure','p_in');
blp.label('内壁随动压力 p');
% Floquet 周期条件（只在本征步开）
pcs = sd.create('pc_s','PeriodicCondition',1);
pcs.selection.set([bWa0 bWaL]);
pcs.set('PeriodicType','Floquet');
pcs.set('kFloquet',{'0';'0';'kz'});
pcs.label('Floquet 周期（管壁）');

%% ---------- acpr ----------
ac = model.physics.create('acpr','PressureAcoustics','geom1');
ac.selection.set(dFl);
pca = ac.create('pc_a','PeriodicCondition',1);
pca.selection.set([bFl0 bFlL]);
pca.set('PeriodicType','Floquet');
pca.set('kFloquet',{'0';'0';'kz'});
pca.label('Floquet 周期（液体）');

%% ---------- 声固耦合 ----------
asb = model.multiphysics.create('asb1','AcousticStructureBoundary','geom1',1);
asb.selection.set(bFSI);

%% ---------- 网格 ----------
me = model.mesh.create('mesh1','geom1');
mp = me.create('map1','Map');
d = mp.create('dz','Distribution');  d.selection.set([ ...
    mphselectbox(model,'geom1',[-e e; -e L+e],'boundary'), ...
    mphselectbox(model,'geom1',[G.r_i-e G.r_i+e; -e L+e],'boundary'), ...
    mphselectbox(model,'geom1',[G.r_o-e G.r_o+e; -e L+e],'boundary')]);
d.set('numelem', opt.nz);
d = mp.create('drf','Distribution'); d.selection.set([bFl0 bFlL]); d.set('numelem', opt.nf);
d = mp.create('drw','Distribution'); d.selection.set([bWa0 bWaL]); d.set('numelem', opt.nw);
me.run;
st = mphmeshstats(model,'mesh1');
say('网格：%d 四边形（液 %d 层 / 壁 %d 层 / 轴向 %d）', st.numelem(strcmp(st.types,'quad')), opt.nf, opt.nw, opt.nz);

%% ---------- Study ----------
model.study.create('std1');
s1 = model.study('std1').create('stat','Stationary');
s1.set('geometricNonlinearity', true);
s1.set('activate', {'solid','on','acpr','off','frame:spatial1','on','frame:material1','on'});
s1.set('useadvanceddisable', true);
s1.set('disabledphysics', {'solid/pc_s'});
s2 = model.study('std1').create('eig','Eigenfrequency');
s2.set('geometricNonlinearity', true);
s2.set('neigs', 25);
s2.set('shift', '1[kHz]');
s2.set('activate', {'solid','on','acpr','on','frame:spatial1','on','frame:material1','on'});
s2.set('useadvanceddisable', true);
s2.set('disabledphysics', {'solid/symz0','solid/blz'});
model.study('std1').setGenPlots(false);  model.study('std1').setGenConv(false);
model.sol.create('sol1');  model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');  model.sol('sol1').createAutoSequence('std1');

say('T4 模型建好：%s（%s，%s）', char(model.tag), opt.wallmodel, opt.fluid);
end
