%% 绘制模式分布图
clear;
clc;
close all;

%% 定义参数
a = 60 * 1e-3; % 横截面宽度a
b = 5.08 * 1e-3; % 横截面宽度b
c = 3e8; % 光速

%% 计算截止频率TEmn
fcmn = zeros(2, 2);

for m = 1:2

    for n = 0:1
        lamdac = 2 * pi / sqrt((m * pi / a) ^ 2 + (n * pi / b) ^ 2);
        fc = c / lamdac;
        fcmn(m, n + 1) = fc;
    end

end

%% 绘图
figure;
hold on

% 确定最大和最小的截止频率
max_fc = max(fcmn(:));
min_fc = min(fcmn(:));

% 计算线段长度的比例系数
length_coeff = 1 / (max_fc - min_fc);

% 遍历所有模式并绘制
for i = 1:2

    for j = 1:2
        fc = fcmn(i, j);
        length = (max_fc - fc) * length_coeff; % 计算线段长度

        % 绘制线段
        plot([0, length], [fc, fc] / 1e9, 'Color', 'b', 'LineWidth', 1);

        % 绘制圆圈
        plot(length, fc / 1e9, 'ro', 'Color', 'b', 'MarkerSize', 4);

        % 添加文字标注
        text(length + 0.05, fc / 1e9, sprintf('TE%d%d', i, j - 1), 'HorizontalAlignment', 'left');
    end

end

% 设置图表属性
ylabel('Cutoff Frequency (GHz)');
title('Mode Distribution Plot');
grid on;
xlim([0 1.2])
ylim([0 10])
hold off
