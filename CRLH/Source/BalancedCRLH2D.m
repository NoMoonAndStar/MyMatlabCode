%% 计算二维的CRLH的色散曲面

%初始化
clear
clc
close all

%四个电路参数
LR = 2.5 * 1e-9; %右手单位长度电感
CR = 1e-12; %右手单位长度电感
LL = 2.5 * 1e-9; %右手单位长度电感
CL = 1e-12; %右手单位长度电感

%扫频范围
Omega = 1e11;
omega = 0:Omega / 100000:Omega;

%晶格间距
px = 10 * 1e-3;
py = 20 * 1e-3;

%电路
Z = 1i * omega * LR / 2 + 1 ./ (1i * 2 * omega * CL);
Y = 1i * omega * CR + 1 ./ (1i * omega * LL);
x = -Z .* Y;

%gamma->X
kpx = acos(1 - x / 2);
kpy = 0;
%这里把虚部为0的点去掉
idx = imag(kpx) == 0;
kpx = kpx(idx);
omega1 = omega(idx);

%因为Beta是向量，所以取范数
Beta = sqrt(real(kpy) .^ 2 + real(kpx) .^ 2);

h1 = figure;
subplot(1, 3, 1)
plot(Beta, omega1, '.', 'Color', 'b')
ylim([0 Omega])
xlabel('Beta*p')
ylabel('omega')
title('Dispersion relation from gamma to X')

%X->M
kpx = pi;
kpy = acos(3 - x / 2);
idx = imag(kpy) == 0;
kpy = kpy(idx);
omega2 = omega(idx);

Beta = sqrt(real(kpy) .^ 2 + real(kpx) .^ 2);

subplot(1, 3, 2)
plot(Beta, omega2, '.', 'Color', 'b')
ylim([0 Omega])
xlabel('Beta*p')
ylabel('omega')
title('Dispersion relation from X to M')

%M->gamma
kpx = acos(1 - x / 4);
kpy = acos(1 - x / 4);
idx = imag(kpx) == 0;
kpx = kpx(idx);
idx = imag(kpy) == 0;
kpy = kpy(idx);

omega3 = omega(idx);

Beta = sqrt(real(kpy) .^ 2 + real(kpx) .^ 2);
%这里反转Beta是为了便于画图，没有实际意义
Beta = -Beta + max(Beta);

subplot(1, 3, 3)
plot(Beta, omega3, '.', 'Color', 'b')
ylim([0 Omega])
xlabel('Beta*p')
ylabel('omega')
title('Dispersion relation from M to gamma')
