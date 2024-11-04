%初始化
clear
clc
close all

%四个电路参数
LR = 2.5 * 1e-9; %右手单位长度电感
CR = 1e-12; %右手单位长度电感
LL = 2.5 * 1e-9; %右手单位长度电感
CL = 1e-12; %右手单位长度电感

Omega = 1e11;
omega = 0:Omega / 20000:Omega; %扫频范围

p = 3.75e-4; %网格间距

Z = 1i * omega * LR / 2 + 1 ./ (1i * 2 * omega * CL);
Y = 1i * omega * CR + 1 ./ (1i * omega * LL);
omega0 = (LR * CR * LL * CL) ^ (-1/4);

gamma = 1 / p * acosh(1 - Z .* Y);
beta = real(gamma);
alpha = imag(gamma);
beta(omega < omega0) = -beta(omega < omega0);

%%fig3.35(b)
h1 = figure;
hold on
plot(alpha * p, omega, '--', 'color', 'r')
plot(beta * p, omega, '.', 'color', 'b')
plot([0 0], [0 omega(end)], 'LineWidth', 1, 'Color', 'k')
xlim([-pi pi])
ylim([0 omega(end)])
xlabel('beta*p or alpha*p')
ylabel('omega')
title('Dispersion relation of the Balanced periodic structure')
legend('alpha', 'beta')
