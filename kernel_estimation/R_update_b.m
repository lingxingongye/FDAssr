function   R    =   R_update_b(MSI_BS, HSI_2D,R,mu)
% 输入：y (向量，大小为 m x 1)，对应的 x (矩阵，大小为 m x n)
% 输出：估计的光谱响应函数 phi_v 的离散值
for i = 1:size(MSI_BS,1)
    y = MSI_BS(i,:)';
    x = HSI_2D';
    [m,n]=size(x);
   

% B 样条基函数的设置
nBasis = n; % B 样条基函数数量
degree = 4; % 
knots = augknt(linspace(1, n, nBasis - degree + 2), degree); % for pavia 

% 构建 B 样条基并评估每个基函数
basisValues = zeros(nBasis, n); % 创建 nBasis x n 的矩阵用于存储基函数值

for j = 1:nBasis
    % 构建每个单独的 B 样条基函数
    Bi = spmak(knots, (1:nBasis == j)); % 生成第 i 个 B 样条基函数
    basisValues(j, :) = fnval(Bi, 1:n);% 评估第 i 个样条基函数在波段点上的值
    
end
% 构建线性系统
% 目标是求解：y = X * c，其中 X 是组合后的 x 和 basisValues
X = zeros(m, nBasis);
for l = 1:m
    X(l, :) = x(l, :) * basisValues';
end

% 构建二次导数平滑惩罚矩阵 P
    dbasis = fnder(spmak(knots, eye(nBasis)), 2); % 对 B 样条基函数求二阶导
    secondDerivValues = fnval(dbasis, 1:n); % 评估二阶导数值
    P = secondDerivValues * secondDerivValues'; % 惩罚矩阵，基于二阶导数的内积

  %  构造目标函数的矩阵 H 和向量 f
    lambda = 0.0005;
    H = X' * X + lambda * P;
    f = -X' * y;
%     % 求解带平滑正则化的优化问题
    options = optimoptions('quadprog', 'Display', 'off');
    c = quadprog(H, f, [], [], [], [], zeros(size(X, 2), 1), [], [], options);

 % 在估计系数 c 时施加约束，使得 phi_v >= 0 且 sum(phi_v) = 1
   % 对 phi_v 施加约束，将 phi_v = basisValues' * c 展开
%     Aeq = sum(basisValues, 2)'; % 等式约束，确保 phi_v 的和为 1
%     beq = 1;
%     A_ineq = -basisValues'; % 确保 phi_v >= 0
%     b_ineq = zeros(n, 1);
% % 
% % %     % 使用 lsqlin 求解带有约束的优化问题
%     options = optimoptions('lsqlin', 'Display', 'off');
%     c = lsqlin(X, y, A_ineq, b_ineq, Aeq, beq, [], [], [], options);
%       c = lsqlin(X, y, A_ineq, b_ineq, [], [], [], [], [], options);
% 求解 B 样条系数 c
% c = lsqnonneg(X, y);
c = lsqlin(X, y);

% 计算估计的光谱响应函数 phi_v
phi_v = basisValues' * c;   

R(i,:) = phi_v;
end
end




