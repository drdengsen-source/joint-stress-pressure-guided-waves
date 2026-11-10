function D = tensor_to_D5(A)
% TENSOR_TO_D5  把 3×3×3×3 张量（基 (r,θ,z) = (1,2,3)）约化为轴对称无扭转的 5×5 矩阵
%   梯度分量顺序：h = [H_rr, H_θθ, H_zz, H_rz, H_zr]，D(a,b) = A(ia,ja,ib,jb)
%   增量应力 δP_(ij) = D h；能量密度 ½ h' D h。D 有大对称（对称矩阵），无小对称。
idx = [1 1; 2 2; 3 3; 1 3; 3 1];
D = zeros(5,5);
for a = 1:5
    for b = 1:5
        D(a,b) = A(idx(a,1), idx(a,2), idx(b,1), idx(b,2));
    end
end
end
