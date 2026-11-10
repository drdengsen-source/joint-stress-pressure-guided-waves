function K = coilK(zc, Nt, sgn, rc, rsurf, pitch, z)
%COILK  周向单层线圈在管外壁感应的面电流密度（良导体镜像近似，每安培线圈电流）
%   K = coilK(zc, Nt, sgn, rc, rsurf, pitch, z)
%   线圈 Nt 匝单层、匝距 pitch、导线中心半径 rc、轴向中心 zc、电流 sgn×1 A（+φ 为正）。
%   K_φ(z) = −2 H_z,coil(r_surf, z)：镜像使切向 H 加倍，面电流 K = n × H。
%   单位 A/m 每 A。适用条件：趋肤深度 ≪ 离壁距离 ≪ 管半径
%   （300–520 kHz 时 7075 铝 δ = 0.21–0.16 mm；离壁 0.5 mm；r_o = 6.5 mm）。
mu0 = 4e-7*pi;
zt  = zc + ((1:Nt) - (Nt+1)/2)*pitch;
K   = 0;
for k = 1:Nt
    [~, bz] = loopB(rc, zt(k), sgn, rsurf, z);
    K = K - 2*bz/mu0;
end
end
