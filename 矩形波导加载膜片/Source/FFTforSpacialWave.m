%% 对空间波形进行傅里叶变换
clear;
clc;
close all;

%% 导入数据
Unloaded = load('..\Data\无加载空间波形.txt');
Loaded = load('..\Data\离散加载空间波形18GHz.txt');

%% Post-processing取一个波长的数据从波峰到波峰
Ts1 = (Unloaded(2, 1) - Unloaded(1, 1)) * 1e-3;
Ts2 = (Loaded(2, 1) - Loaded(1, 1)) * 1e-3;
Fs1 = 1 / Ts1;
Fs2 = 1 / Ts2;
f1 = Unloaded(:, 2);
f2 = Loaded(:, 2);
peak1 = [1007, 2450];
% peak2 = [761, 1881];  % 3GHz离散
% peak2 = [1541, 5422]; % 3GHz密集
% peak2 = [55, 855];    % 8GHz离散
peak2 = [261, 941]; % 18GHz离散
t01 = (0:(peak1(2) - peak1(1)))' * Ts1;
t02 = (0:(peak2(2) - peak2(1)))' * Ts2;
f1 = f1(peak1(1):peak1(2));
f2 = f2(peak2(1):peak2(2));

%% 对波形进行重复后傅里叶变换
N = 1000;
f1 = repmat(f1, N, 1);
f2 = repmat(f2, N, 1);
t1 = t01;
t2 = t02;

for i = 2:N
    t01 = t01 + (peak1(2) - peak1(1)) * Ts1;
    t02 = t02 + (peak2(2) - peak2(1)) * Ts2;
    t1 = [t1; t01];
    t2 = [t2; t02];
end

L1 = length(f1);
L2 = length(f2);
fx1 = Fs1 / L1 * (0:(L1 / 2));
fx2 = Fs2 / L2 * (0:(L2 / 2));
y1 = fft(f1);
y2 = fft(f2);

%% 绘图
h1 = figure;
hold on
plot(t1, f1, 'lineWidth', 3')
plot(t2, f2, 'lineWidth', 3')
xlabel('Distance(m)')
ylabel('Real(E)(V/m)')
legend('无加载', '有加载')
title('Space Waveform')

h2 = figure;
P21 = abs(y1 / L1);
P22 = abs(y2 / L2);
P11 = P21(1:L1 / 2 + 1);
P12 = P22(1:L2 / 2 + 1);
P11(2:end - 1) = 2 * P11(2:end - 1);
P12(2:end - 1) = 2 * P12(2:end - 1);
%转每米为弧度
fx1 = fx1 * 2 * pi;
fx2 = fx2 * 2 * pi;
%归一化
P11 = P11 / sum(P11);
P12 = P12 / sum(P12);
hold on
plot(fx1, P11, 'lineWidth', 3')
plot(fx2, P12, 'lineWidth', 3')
xlim([0 4000])
xlabel('Beta(rad/m)')
ylabel('magnitude')
legend('无加载', '有加载')
title('FFT of Space Waveform(Normalized)')
