function resultMatrix = Remap(omega, k)

    % 初始化一个空单元格数组，用于存储每个扩展后的块矩阵
    blocks = {};
    
    % 遍历向量omega中的每个元素
    for i = 1:length(omega)
        % 生成k x k的矩阵，所有元素都等于omega(i)
        block = ones(k, k) * omega(i);
        
        % 将生成的块矩阵添加到单元格数组中
        blocks{i} = block;
    end
    
    % 使用blkdiag函数将所有块矩阵沿对角线堆叠成一个大矩阵
    resultMatrix = blkdiag(blocks{:});

end