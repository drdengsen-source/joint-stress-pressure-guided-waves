function c = lamb_S0_phase(fd, cL, cT)
% LAMB_S0_PHASE  无预应力各向同性板 S0 模态相速度（Rayleigh–Lamb 精确解）
%   fd：频厚积 [Hz·m]；cL, cT 体波速。返回 c [m/s]。
%   对称模态特征方程（cT < c < cL，p^2 < 0 写成双曲形式）：
%   (q^2 − k^2)^2 cosh(αh) sin(qh) − 4 k^2 α q sinh(αh) cos(qh) = 0
%   α = k sqrt(1 − c^2/cL^2),  q = k sqrt(c^2/cT^2 − 1),  h = d/2
d = 1e-3;  f = fd/d;  h = d/2;  w = 2*pi*f;
F = @(c) charS0(c, w, h, cL, cT);
c = fzero(F, [cT*1.001, cL*0.999]);
end

function y = charS0(c, w, h, cL, cT)
k = w/c;
al = k*sqrt(1 - c^2/cL^2);
q  = k*sqrt(c^2/cT^2 - 1);
y = (q^2 - k^2)^2*cosh(al*h)*sin(q*h) - 4*k^2*al*q*sinh(al*h)*cos(q*h);
% 归一化，避免 k^4 量级过大
y = y/k^4;
end
