function dk = safe_perturb(S, md, omega, dD, dKs, drho_rr, drho_zz, drho_s, Ks, rho_f)
% SAFE_PERTURB  T2 §2.1 的一阶摄动：Δk|_ω = (ω/P)(Δ⟨T⟩ − Δ⟨V⟩)_场固定
%   dD      : 5×5 或 5×5×nGs，壁增量模量 ΔA（位移梯度形式）
%   dKs     : 液体绝热体积模量增量（含拉回几何项）
%   drho_rr, drho_zz : 液体（各向异性）密度增量，用于 Δ⟨T⟩_f = ¼∫(Δρ_rr|v_r|²+Δρ_zz|v_z|²)
%   drho_s  : 壁密度增量（温度用）
%   Ks, rho_f : 未扰动液体参数（算 |p|²/K² 与 |v|² 用）
% 返回 Δk（每单位扰动，1/m）；Δc/c = −Δk/k。
uniformD = (ndims(dD) == 2);
gs = S.s;  dVs = 0;  dTs = 0;
for i = 1:numel(gs.r)
    w = gs.w(i)*gs.r(i);
    if uniformD, dDi = dD; else, dDi = dD(:,:,i); end
    hi = md.h(:,i);  ui = md.u(:,i);
    dVs = dVs + 0.25*w*real(hi'*dDi*hi);
    dTs = dTs + 0.25*w*drho_s*omega^2*(ui'*ui);
end
gf = S.f;  dVf = 0;  dTf = 0;
for i = 1:numel(gf.r)
    w = gf.w(i)*gf.r(i);
    dVf = dVf + 0.25*w*(dKs/Ks^2)*abs(md.p(i))^2;
    dTf = dTf + 0.25*w*(drho_rr*abs(md.dp(i))^2 + drho_zz*md.k^2*abs(md.p(i))^2)/(omega^2*rho_f^2);
end
dk = (omega/md.Pw)*(dTs + dTf - dVs - dVf);
end
