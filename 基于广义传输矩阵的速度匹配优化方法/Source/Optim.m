%% 优化电路————终极无敌版
clear
clc
close all

%% 定义参数区域！
UnitE = 1 * pi / 180;
F = 1e9;
VaractorCj0 = 8.47e-12;
VaractorM = 70;
VaractorVj0 = 80;
VStatic = 3;
Zc1 = 50;
Zc2 = 50;
LSeries1 = 1e-9;
LSeries2 = 1e-9;
ZPort1 = 50 + 1i * 50;
ZPort2 = 50 + 1i * 50;
ZPort3 = 50;
ZPort4 = 50;
ZPort5 = 50;
ZPort6 = 50;
N = 24;
Fmin = 0;
Fmax = 1e9;
FSweep = 1001;
Freq = linspace(Fmin, Fmax, FSweep); % 扫频范围

%% 计算S参数
tic
S = CaculateSParameter(UnitE, F, VaractorCj0, VaractorM, VaractorVj0, VStatic, Zc1, Zc2, ...
    LSeries1, LSeries2, ZPort1, ZPort2, ZPort3, ZPort4, ZPort5, ZPort6, N, Fmin, Fmax, FSweep);
toc

% 提取S41,S52,S44,S55
S41 = squeeze(S(4, 1, :));
S52 = squeeze(S(5, 2, :));
S44 = squeeze(S(4, 4, :));
S55 = squeeze(S(5, 5, :));

dBS41 = 20 * log10(abs(S41));
dBS52 = 20 * log10(abs(S52));
dBS44 = 20 * log10(abs(S44));
dBS55 = 20 * log10(abs(S55));

% 导入仿真S参数
SimSParameter = importdata('..\Data\SParameter_Sim.csv');
SimdBS41 = SimSParameter(:, 3);
SimdBS52 = SimSParameter(:, 4);
SimdBS44 = SimSParameter(:, 1);
SimdBS55 = SimSParameter(:, 2);

% 绘图，与仿真对比是否正确
figure(1);
subplot(2, 2, 1)
hold on
plot(Freq, dBS41)
plot(Freq, SimdBS41)
legend('Theory', 'Sim')
xlabel('Frequency(Hz)')
ylabel('dB')
title('|S41|')
grid on
subplot(2, 2, 2)
hold on
plot(Freq, dBS52)
plot(Freq, SimdBS52)
legend('Theory', 'Sim')
xlabel('Frequency(Hz)')
ylabel('dB')
title('|S52|')
grid on
subplot(2, 2, 3)
hold on
plot(Freq, dBS44)
plot(Freq, SimdBS44)
legend('Theory', 'Sim')
xlabel('Frequency(Hz)')
ylabel('dB')
title('|S44|')
grid on
subplot(2, 2, 4)
hold on
plot(Freq, dBS55)
plot(Freq, SimdBS55)
legend('Theory', 'Sim')
xlabel('Frequency(Hz)')
ylabel('dB')
title('|S55|')
grid on
