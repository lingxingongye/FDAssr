function   B    =   B_update_b2(R_HSI_up, MSI,size_B,sf,B,mu)
aaa = [];
bbb = [];

for i= 1*sf+1:sf:size(MSI,1)- floor((size_B-1)/2)-1 
    for j=1*sf+1:sf:size(MSI,2)- floor((size_B-1)/2)-1
        for k=1:size(MSI,3)
            d_size=floor((size_B-1)/2);
      image_patch=MSI(i-d_size:i+d_size,j-d_size:j+d_size,k);
       yyy =  image_patch(:); 
        E =  R_HSI_up(i, j, k);
        aaa = [aaa;yyy'];
        bbb = [bbb;E];
        end
    end   
end
[m,n]=size(aaa);
% B 样条基函数的设置
nBasis = n; % B 样条基函数数量
degree = 3; % B 样条的阶次（即三次 B 样条）
knots = augknt(linspace(1, n, nBasis - 1), degree); % for pavia 
% knots = augknt([1 1 linspace(1, n, nBasis - degree) n n], degree); %for Houston

% 构建 B 样条基并评估每个基函数
basisValues = zeros(nBasis, n); % 创建 nBasis x n 的矩阵用于存储基函数值
for j = 1:nBasis
    % 构建每个单独的 B 样条基函数
    Bi = spmak(knots, (1:nBasis == j)); % 生成第 i 个 B 样条基函数
    basisValues(j, :) = fnval(Bi, 1:n); % 评估第 i 个样条基函数在波段点上的值
end

% 构建线性系统
% 目标是求解：y = X * c，其中 X 是组合后的 x 和 basisValues
X = zeros(m, nBasis);
for l = 1:m
    % x(i, :) 是 1 x n 的行向量，basisValues 是 nBasis x n，
    % 所以 x(i, :) * basisValues' 的大小为 1 x nBasis
    X(l, :) = aaa(l, :) * basisValues';
end

% 求解 B 样条系数 c
c = lsqnonneg(X, bbb);
% c = lsqlin(X, bbb);

phi_v = basisValues' * c;


B=reshape(phi_v,[size(B,1),size(B,2)]);

end



