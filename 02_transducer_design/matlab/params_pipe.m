function P = params_pipe(grade, fluid)
% PARAMS_PIPE  充液管道理论计算的参数表（换牌号／换液体只改这里）
%   P = params_pipe            默认 '7075' + 'oil'
%   P = params_pipe('7075','water')
% 单位：SI（m, kg, s, Pa, K）。应力拉为正。
%
% 参数取值与出处：
%   管材    7075-T6 无缝管，外径 13 mm，壁厚 1 mm
%   液体    矿物基液压油（ISO VG 46），设计值 c_f = 1450 m/s、ρ_f = 865 kg/m³ @20 °C
%   B/A     取 10（烃类液体文献范围 9–11）
% 出处：
%   几何    外径 13 mm，壁厚 1 mm
%   7075    λ, μ, l, m, n：Stobbe D.M. (2005) M.S. thesis, Georgia Tech,
%           “Acoustoelasticity in 7075-T651 aluminum and dependence of third
%           order elastic constants on fatigue damage”。数值转引自
%           Gandhi, Michaels & Lee, J. Acoust. Soc. Am. 132 (2012) 1284, Table I。
%           ρ = 2810 为 7075 手册值。
%   6061    仅列二阶常数，本文未使用。
%   液体    ρ_f：Castrol Hyspin VG46 数据表 0.870 g/mL @15 °C；Heliyon 8 (2022) e11831
%           实测 HL32 油 865.5 kg/m³ @20 °C，ρ = ρ15[1 − 0.0007(T − 15)] → β_f = 7.0e-4/K。
%           c_f：矿物液压油纯油约 1400–1500 m/s（Hawe 流体手册"纯油≈1500，含气取 1300"），
%           设计值 1450；实测值见正文 4.1 节。
%           dc_f/dT：烃类液体典型 −3 ~ −4 m/s/K，取 −3.5（实测 −3.52，见补充材料表 S8）。
%           B/A：烃类/油脂类液体文献范围 9–11（Beyer, Nonlinear Acoustics 汇编；生物脂类 10–11），
%           取 10。水：B/A = 5.0 @20 °C（Beyer 1960）。
%           c_p ≈ 1900 J/(kg K)、γ = K_s/K_T ≈ 1.10：矿物油典型量级。
%   温度    (1/E)dE/dT、α 为铝合金典型量级；温度系数的实测值见正文 4.2 节。

if nargin < 1 || isempty(grade), grade = '7075'; end
if nargin < 2 || isempty(fluid), fluid = 'oil'; end

%% 几何
G.r_o = 6.5e-3;  G.r_i = 5.5e-3;  G.t = 1.0e-3;  G.r_m = 6.0e-3;
G.W_E = 0.100;                                  % EMAT 标距（空间固定）
G.A_wall = pi*(G.r_o^2 - G.r_i^2);              % 37.70 mm^2
G.A_bore = pi*G.r_i^2;                          % 95.03 mm^2
G.sig_theta_per_p = G.r_i/G.t;                  % 厚度平均环向应力 / p  (= 5.5)
G.sig_z_closed_per_p = G.A_bore/G.A_wall;       % 闭端轴向应力 / p      (= 2.52)
G.sig_r_per_p = -0.5;                           % 径向应力厚度平均 / p（Lamé 线性分布的均值，估计用）
P.geom = G;

%% 管壁材料表
switch upper(grade)
    case '7075'
        W.name = 'Al 7075-T651 (Stobbe 2005, via Gandhi 2012)';
        W.rho  = 2810;              % kg/m^3，7075 手册值
        W.lam  = 54.308e9;  W.mu = 27.174e9;
        W.l    = -252.2e9;  W.m  = -324.9e9;  W.n = -351.2e9;
    case '6061'
        W.name = 'Al 6061-T6 (second-order constants only)';
        W.rho  = 2700;
        W.lam  = 51.1e9;    W.mu = 26.3e9;      % E≈69.6 GPa, ν≈0.33
        W.l    = NaN;       W.m  = NaN;         W.n = NaN;
    otherwise
        error('params_pipe: 未知牌号 %s', grade);
end
W.E   = W.mu*(3*W.lam + 2*W.mu)/(W.lam + W.mu);
W.nu  = W.lam/(2*(W.lam + W.mu));
W.K0  = W.lam + 2*W.mu/3;                       % 体积模量
W.cL  = sqrt((W.lam + 2*W.mu)/W.rho);
W.cT  = sqrt(W.mu/W.rho);
W.c_plate = sqrt(W.E/((1 - W.nu^2)*W.rho));     % 薄板 S0 极限速度
W.c_bar   = sqrt(W.E/W.rho);                    % 杆速（Korteweg 极限下的环频用它）
W.alpha   = 23.6e-6;                            % 线膨胀系数 1/K（典型值）
W.dlnE_dT = -5.0e-4;                            % (1/E)dE/dT, 1/K（典型量级）
W.dnu_dT  = 0;                                  % 取 0
P.wall = W;

%% 液体表
switch lower(fluid)
    case 'oil'
        F.name  = '矿物液压油（设计值）';
        F.c     = 1450;     F.rho = 865;    F.BA = 10;
        F.beta  = 7.0e-4;                   % 体膨胀系数 1/K
        F.cp    = 1900;                     % J/(kg K)
        F.dc_dT = -3.5;                     % m/s/K（矿物油典型值）
        F.gamma = 1.10;                     % c_p/c_v；K_T = K_s/γ
    case 'water'
        F.name  = '水 20 °C（文献典型值）';
        F.c     = 1482;     F.rho = 998;    F.BA = 5.0;
        F.beta  = 2.1e-4;   F.cp = 4182;    F.dc_dT = +2.4;   F.gamma = 1.007;
    otherwise
        error('params_pipe: 未知液体 %s', fluid);
end
F.T0  = 293.15;
F.Ks  = F.rho*F.c^2;                            % 绝热体积模量（声波用）
F.KT  = F.Ks/F.gamma;                           % 等温体积模量（静态加压用）
F.a_T = F.dc_dT/F.c;                            % (1/c)(∂c/∂T)_p
% 压力系数 a_p = (1/c)(∂c/∂p)：
%   绝热：a_p = (B/2A)/K_s
%   等温（Beyer 关系 B/A = 2ρc(∂c/∂p)_T + (2cTβ/c_p)(∂c/∂T)_p）：
F.a_p_s = (F.BA/2)/F.Ks;
F.a_p   = (F.BA - 2*F.c^2*F.T0*F.beta*F.a_T/F.cp)/(2*F.Ks);
P.fluid = F;

P.grade = grade;  P.fluidName = fluid;
end
