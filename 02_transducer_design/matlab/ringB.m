function [Br, Bz] = ringB(type, mg, zc, sgn, M, r, z)
%RINGB  均匀磁化环形永磁体的静磁场（等效面电流叠加 loopB）
%   [Br, Bz] = ringB(type, mg, zc, sgn, M, r, z)
%   type = 'R' 径向充磁（sgn = +1 磁化向外，此时内孔面是 S 极、外柱面是 N 极）：
%            等效电流在两端面，顶面(z = zc + W/2) −M φ̂、底面 +M φ̂，径向电流密度 M (A/m)
%   type = 'A' 轴向充磁（sgn = +1 磁化沿 +z，+z 端面是 N 极）：
%            等效电流在内外柱面，外柱面 +M φ̂、内柱面 −M φ̂（沿 z 的电流密度 M）
%   mg.ri, mg.ro, mg.W：内半径、外半径、轴向宽 (m)；zc 环中心；M = B_rem/μ0 (A/m)；场点 (r, z)。
%   适用：周围介质 μ_r = 1（铝、树脂、空气；NdFeB 回复导磁率 1.05 取 1），磁体外任意点精确。
%   输出 T。离散：径向环每 0.05 mm 一圈，轴向环每 0.1 mm 一圈（离场点 ≥ 0.5 mm 时误差 < 1%）。
Br = 0; Bz = 0;
switch upper(type)
    case 'R'
        na = max(40, round((mg.ro - mg.ri)/0.05e-3));
        ae = linspace(mg.ri, mg.ro, na+1);
        a  = (ae(1:end-1) + ae(2:end))/2;
        da = (mg.ro - mg.ri)/na;
        for k = 1:na
            [b1r, b1z] = loopB(a(k), zc - mg.W/2, +sgn*M*da, r, z);
            [b2r, b2z] = loopB(a(k), zc + mg.W/2, -sgn*M*da, r, z);
            Br = Br + b1r + b2r;   Bz = Bz + b1z + b2z;
        end
    case 'A'
        nz = max(80, round(mg.W/0.1e-3));
        ze = linspace(zc - mg.W/2, zc + mg.W/2, nz+1);
        zs = (ze(1:end-1) + ze(2:end))/2;
        dz = mg.W/nz;
        for k = 1:nz
            [b1r, b1z] = loopB(mg.ro, zs(k), +sgn*M*dz, r, z);
            [b2r, b2z] = loopB(mg.ri, zs(k), -sgn*M*dz, r, z);
            Br = Br + b1r + b2r;   Bz = Bz + b1z + b2z;
        end
    otherwise
        error('ringB: type 只能是 ''R''（径向充磁）或 ''A''（轴向充磁）');
end
end
