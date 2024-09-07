CSTData = load('..\Data\色散曲线第二种.txt');
TheoryData = load('..\Data\Disepersion.mat');

h1 = figure;
hold on
%从CST仿真绘制色散曲线
plot(CSTData(:, 1), CSTData(:, 2), 'bo', 'MarkerSize', 4)

%从理论计算绘制色散曲线
plot(TheoryData.Beta(:, 1) * 180 / pi, TheoryData.k, 'r.', 'LineWidth', 1)

legend('CST', 'Theory')
ylim([0 25])
xlabel('Beta*p(deg)')
ylabel('Frequency(GHz)')
title('The dispersion curve of the periodic structure')
