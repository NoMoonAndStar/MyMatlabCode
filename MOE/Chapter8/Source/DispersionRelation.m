%%8.1节色散关系
d = 0.01;
syms betad

%真空中的波数
kd = 0:0.001:4;
%加载后的相位常数
realbetad = zeros(1, length(kd));

for i = 1:length(kd)
    %计算相位常数
    f = cos(betad) - cos(kd(i)) + 2 * kd(i) * sin(kd(i));
    s = solve(f == 0);

    %判断解是否为实数
    if isreal(s)

        if sum(double(s) >= 0)
            realbetad(i) = min(double(s(double(s) >= 0)));
        end

    end

end

%剔除虚数解，除第一个0外都删除
realbetad2 = realbetad(2:end);
kd2 = kd(2:end);
idxreal = realbetad2 > 0;
realbetad = [0, realbetad2(idxreal)];
kd = [0, kd2(idxreal)];

%绘制一边的色散曲线
plot(realbetad, kd, '.');
xlabel('betad')
ylabel('kd')
title('Dispersion relation')

%绘制另一半的色散曲线
hold on
plot(-realbetad, kd, '.');
