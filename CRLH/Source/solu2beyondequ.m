function [x_solutions,y_solutions] = solu2beyondequ(a,b,c,m)

    kpx = linspace(-pi, pi, m);
    kpy = linspace(-pi, pi, m);
    
    % 初始化一个空数组来存储解
    x_solutions = [];
    y_solutions = [];
    
    % 遍历x和y的范围
    for x = kpx
        for y = kpy
            % 计算方程值
            equation_value = imag(a*cos(x) + b*cos(y) + c);
            
            % 如果方程值接近于0（考虑浮点数精度），则认为找到一个解
            if abs(equation_value) < 5e-2
                x_solutions = [x_solutions; x];
                y_solutions = [y_solutions; y];
            end
        end
    end
end
