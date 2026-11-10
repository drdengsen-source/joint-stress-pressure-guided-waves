function E4_paperfig_data()
%E4_PAPERFIG_DATA  给论文图 A 导出数据：D1T 叠列在管外壁的 B_r、B_z（E1）与两组线圈最终设计的力分布（E2）
% 输出 E4_paperfig_data.csv：z_mm, Br_T, Bz_T, fz_wall_N_per_m_A, fr_liq_N_per_m_A, fr_wall_x, fz_liq_x
here = fileparts(mfilename('fullpath'));
F = load(fullfile(here, 'E1_fields.mat'), 'z', 'D1Tw');
C = load(fullfile(here, 'E2_coils.mat'), 'sweep');
S = C.sweep.St{3};                                   % D1T
z = F.z(:);
M = [z*1e3, F.D1Tw.Br(:), F.D1Tw.Bz(:), S.fzw(:), S.frl(:), S.frw(:), S.fzl(:)];
fid = fopen(fullfile(here, 'E4_paperfig_data.csv'), 'w');
fprintf(fid, 'z_mm,Br_T,Bz_T,fz_wall,fr_liq,fr_wall_x,fz_liq_x\n');
fprintf(fid, '%.3f,%.6f,%.6f,%.4f,%.4f,%.4f,%.4f\n', M.');
fclose(fid);
fprintf('exported %d rows; wall coil %.2f mm x %d turns, liquid coil %.2f mm x %d turns\n', size(M, 1), ...
    C.sweep.ww(1)*0 + 5.25, 35, 2.85, 19);
end
