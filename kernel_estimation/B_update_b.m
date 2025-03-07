function   B    =   B_update_b(R_HSI_up, MSI,size_B,sf,B,mu)

  X = MSI;
  Y = R_HSI_up;

    [M, N,~] = size(Y); % 获取输入图像 Y 的尺寸
    kernel_size = size_B; % 卷积核的大小为 7x7
    d_size=floor((size_B-1)/2);% 计算卷积核的一半大小，用于邻域窗口提取
    % 构建B样条基函数
    % 构建水平和竖直方向的基函数
    nBasis = kernel_size; % B样条基函数数量，设为 kernel_size - 1，即 6 个基函数
    degree = 3;
%         nBasis = kernel_size; % B样条基函数数量，设为 kernel_size - 1，即 6 个基函数
%     degree = 4;
    knots = augknt(linspace(0, 1, nBasis - degree + 2), degree); % 构建节点矢量，节点用于构造B样条基函数

    
    Bx = zeros(kernel_size, nBasis); % 初始化水平方向的B样条基函数矩阵，大小为 7 x (nBasis)
    By = zeros(kernel_size, nBasis); % 初始化竖直方向的B样条基函数矩阵，大小为 7 x (nBasis)
    
    for i = 1:nBasis
        % 构建每个单独的B样条基函数
        Bi_x = spmak(knots, (1:nBasis == i)); % 构建第 i 个水平方向的B样条基函数
        Bi_y = spmak(knots, (1:nBasis == i)); % 构建第 i 个竖直方向的B样条基函数
        Bx(:, i) = fnval(Bi_x, linspace(0, 1, kernel_size)); % 评估水平方向的B样条基函数值，存储在 Bx 中
        By(:, i) = fnval(Bi_y, linspace(0, 1, kernel_size)); % 评估竖直方向的B样条基函数值，存储在 By 中
    end
    
    % 构造矩阵 A 和 向量化后的 X
    A = []; % 初始化矩阵 A，用于存储线性系统的系数矩阵
    X_vector = []; % 初始化向量 X_vector，用于存储目标输出值
     
    

    for m = 1*sf+1:sf:size(X,1)- floor((size_B-1)/2)-1
        for n = 1*sf+1:sf:size(X,2)- floor((size_B-1)/2)-1
            for k=1:size(X,3)  
            % 取 Y 的邻域窗口
            patch = X(m-d_size:m+d_size, n-d_size:n+d_size,k); % 提取当前像素点的 7x7 邻域窗口
            % 使用 B 样条基函数展开
            A_row = zeros(nBasis * nBasis, 1);
            idx = 1;
            for i = 1:nBasis
                for j = 1:nBasis
                        A_row(idx) = sum(sum(patch .* (Bx(:, i) * By(:, j)')));
                        idx = idx + 1;
                end
            end
            A = [A; A_row']; % 将当前行添加到矩阵 A 中
            X_vector = [X_vector; Y(m, n, k)]; % 将对应位置的输出值添加到向量 X_vector 中
            end
        end
    end
  
    
% 构造线性等式约束 A_eq_B 和 b_eq，用于约束最终 B 的元素和为 1
A_eq_B = zeros(1, nBasis * nBasis); % 初始化 A_eq_B
for i = 1:kernel_size
    for j = 1:kernel_size
        % 对每个卷积核元素 B(i, j)，找出与之对应的基函数系数组合
        idx = 1;
        for p = 1:nBasis
            for q = 1:nBasis
                A_eq_B(idx) = A_eq_B(idx) + (Bx(i, p) * By(j, q));
                idx = idx + 1;
            end
        end
    end
end

% 约束条件为 B 的元素和为 1
b_eq = 1;

% 构造线性不等式约束 A_ineq 和 b_ineq，使得 c_coeff 的所有元素 >= 0
A_ineq = -eye(nBasis * nBasis); % 每个系数 >= 0
b_ineq = zeros(nBasis * nBasis, 1); % 右侧全为 0，表示非负约束

% 使用 lsqlin 求解
options = optimoptions('lsqlin', 'Display', 'off'); % 关闭显示
c_coeff = lsqlin(A, X_vector, A_ineq, b_ineq, A_eq_B, b_eq, [], [], [], options);


    c_coeff_matrix = reshape(c_coeff, [nBasis, nBasis]);
    % 计算估计的卷积核
    c = zeros(kernel_size, kernel_size); % 初始化卷积核矩阵 c，大小为 7x7
    for i = 1:kernel_size
        for j = 1:kernel_size
           % 使用基函数展开计算卷积核的最终离散值
            c(i, j) = sum(sum(c_coeff_matrix .* (Bx(i, :)' * By(j, :))));
        end
    end

B = c;

end