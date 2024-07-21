d = 0.01;
syms betad
kd = 0:0.01:8;
realbetad = zeros(1, length(kd));

for i = 1:length(kd)
    f = cos(betad) - cos(kd(i)) + 2 * kd(i) * sin(kd(i));
    s = solve(f == 0);

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

%绘制色散图
plot(realbetad, kd, '.');
