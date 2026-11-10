function [Br, Bz] = loopB(a, z0, I, r, z)
%LOOPB  圆电流环的静磁场（椭圆积分解析式）
%   [Br, Bz] = loopB(a, z0, I, r, z)
%   环半径 a (m)，环所在平面 z = z0，电流 I (A，+φ 方向为正)；
%   场点 (r, z)：同形数组，或标量与数组的组合（按 MATLAB 广播）。
%   出处：Simpson, Lane, Immer, Youngquist, "Simple analytic expressions for the
%   magnetic field of a circular current loop", NASA/TM-2001-210953（柱坐标 K、E 式）。
%   校核：r→0 时退化为 μ0 I a² / [2(a²+z²)^{3/2}]（E1_force_field.m 里自检）。
%   单位：SI，输出 T。
mu0 = 4e-7*pi;
r   = max(r, 1e-9);                       % 轴上奇点
zz  = z - z0;
s   = a.^2 + r.^2 + zz.^2;
al2 = s - 2*a.*r;                         % α²
be2 = s + 2*a.*r;                         % β²
be  = sqrt(be2);
m   = 1 - al2./be2;                       % k²
[K, E] = ellipke(m);
C   = mu0*I/pi;
Br  = C.*zz./(2*al2.*be.*r).*(s.*E - al2.*K);
Bz  = C./(2*al2.*be).*((a.^2 - r.^2 - zz.^2).*E + al2.*K);
end
