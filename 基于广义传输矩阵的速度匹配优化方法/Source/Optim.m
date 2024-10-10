%% 优化电路————终极无敌版
clear
clc
close all

%% 定义不变参数
UnitE = 1 * pi / 180;
F = 1e9;
VaractorCj0 = 8.47e-12;
VaractorM = 70;
VaractorVj0 = 80;
VStatic = 3;
Zc1 = 70;
Zc2 = 50;
LSeries1 = 1e-9;
LSeries2 = 1e-9;
EHorizontalPort1 = 0;
EVerticalPort1 = 0;
EHorizontalPort2 = 0;
EVerticalPort2 = 0;
EHorizontalPort4 = 0;
EVerticalPort4 = 0;
EHorizontalPort5 = 0;
EVerticalPort5 = 0;
N = 24;
Fmin = 0;
Fmax = 1e9;
FSweep = 1001;
Freq = linspace(Fmin, Fmax, FSweep); % 扫频范围

%% 计算S参数
tic
S = CaculateSParameter(UnitE, F, VaractorCj0, VaractorM, VaractorVj0, VStatic, Zc1, Zc2, ...
    LSeries1, LSeries2, EHorizontalPort1, EVerticalPort1, EHorizontalPort2, EVerticalPort2, ...
    EHorizontalPort4, EVerticalPort4, EHorizontalPort5, EVerticalPort5, N, Fmin, Fmax, FSweep);
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
legend('Theory', 'Sim')
xlabel('Frequency(Hz)')
ylabel('dB')
title('|S41|')
grid on
subplot(2, 2, 2)
hold on
plot(Freq, dBS52)
legend('Theory', 'Sim')
xlabel('Frequency(Hz)')
ylabel('dB')
title('|S52|')
grid on
subplot(2, 2, 3)
hold on
plot(Freq, dBS44)
legend('Theory', 'Sim')
xlabel('Frequency(Hz)')
ylabel('dB')
title('|S44|')
grid on
subplot(2, 2, 4)
hold on
plot(Freq, dBS55)
legend('Theory', 'Sim')
xlabel('Frequency(Hz)')
ylabel('dB')
title('|S55|')
grid on
