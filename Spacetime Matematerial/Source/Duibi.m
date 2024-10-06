optimed = load('../Data/优化后放大倍数.txt');
unoptimed = load('../Data/未优化放大倍数.txt');
t = optimed(:, 1);
youhua = optimed(:, 2);
unyouhua = unoptimed(:, 2);
figure;
hold on
plot(t, youhua, 'r', 'linewidth', 1.5);
plot(t, unyouhua, 'b', 'linewidth', 1.5);
legend('优化后放大倍数', '未优化放大倍数');
xlabel('归一化时间 t/T');
ylabel('放大倍数');
title('放大倍数随光速超材料等效电路中经历时间变化');
grid on
