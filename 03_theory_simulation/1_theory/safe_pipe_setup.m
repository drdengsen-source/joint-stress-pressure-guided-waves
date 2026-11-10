function S = safe_pipe_setup(P, nel_s, nel_f)
% SAFE_PIPE_SETUP  轴对称充液管 SAFE 的径向网格与形函数（T3 用）
%   S = safe_pipe_setup(P, nel_s, nel_f)
%   壁 r_i..r_o 分 nel_s 个二次 Lagrange 单元（每单元 3 节点，2 自由度 u_r,u_z），
%   液 0..r_i 分 nel_f 个二次单元（1 自由度 p）。4 点 Gauss。
%   自由度排序：[u_r1 u_z1 u_r2 u_z2 ... | p1 ... pNf]，液节点从 r=0 到 r_i。
if nargin < 2 || isempty(nel_s), nel_s = 8;  end
if nargin < 3 || isempty(nel_f), nel_f = 40; end
G = P.geom;
rs = linspace(G.r_i, G.r_o, 2*nel_s + 1);
rf = linspace(0,     G.r_i, 2*nel_f + 1);
gp = [-0.861136311594053, -0.339981043584856, 0.339981043584856, 0.861136311594053];
gw = [ 0.347854845137454,  0.652145154862546, 0.652145154862546, 0.347854845137454];
Nfun  = @(x) [x.*(x-1)/2, (1-x).*(1+x), x.*(x+1)/2];
dNfun = @(x) [x-0.5, -2*x, x+0.5];

S.geom = G;  S.rs = rs;  S.rf = rf;  S.nel_s = nel_s;  S.nel_f = nel_f;
S.Ns = numel(rs);  S.Nf = numel(rf);  S.ndof = 2*S.Ns + S.Nf;
S.s = gpdata(rs, nel_s, gp, gw, Nfun, dNfun);
S.f = gpdata(rf, nel_f, gp, gw, Nfun, dNfun);
% 耦合节点：壁第 1 节点 (r_i) 的 u_r 自由度，液最后节点 (r_i) 的 p 自由度
S.i_ur_i = 1;
S.i_p_i  = 2*S.Ns + S.Nf;
end

function g = gpdata(rn, nel, gp, gw, Nfun, dNfun)
ng = numel(gp);  n = nel*ng;
g.r = zeros(n,1); g.w = zeros(n,1); g.N = zeros(n,3); g.dN = zeros(n,3); g.nodes = zeros(n,3);
for e = 1:nel
    nodes = (2*e-1):(2*e+1);
    re = rn(nodes);
    for k = 1:ng
        i = (e-1)*ng + k;
        x = gp(k);  Nx = Nfun(x);  dNx = dNfun(x);
        J = dNx*re(:);
        g.r(i) = Nx*re(:);
        g.w(i) = gw(k)*J;          % ∫ f dr ≈ Σ w f
        g.N(i,:) = Nx;  g.dN(i,:) = dNx/J;
        g.nodes(i,:) = nodes;
    end
end
end
