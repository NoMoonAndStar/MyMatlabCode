%% 复现JOSB的图
clear
clc
close all

%% 参数定义
delta = 0; %偏差
alpha = 0.05; %调制强度
Omega = 0.07 * 1e9; %调制频率
omega = 1 * 1e9; %输入频率
eps1 = 1.45;
eps0 = 8.854e-12;
c = 3e8;
cg = c / sqrt(eps1);
g = Omega / cg; %调制空间频率

%% 计算
X = linspace(0, 2 * pi / g, 5000);
%图A
tau = 37.5;
UA = exp(-g * (2 + delta) * alpha * sin(g * X) * tau - (delta + alpha * cos(g * X)) .^ 2/2 / alpha ^ 2 * (-1 - 2 * alpha * g * tau + exp(2 * alpha * g * tau)) ...
    +delta * (delta + alpha * cos(g * X)) * 2 / alpha ^ 2 * (-1 - alpha * g * tau + exp(alpha * g * tau)));

%图B
tau = 150;
UB = exp(-g * (2 + delta) * alpha * sin(g * X) * tau - (delta + alpha * cos(g * X)) .^ 2/2 / alpha ^ 2 * (-1 - 2 * alpha * g * tau + exp(2 * alpha * g * tau)) ...
    +delta * (delta + alpha * cos(g * X)) * 2 / alpha ^ 2 * (-1 - alpha * g * tau + exp(alpha * g * tau)));

%图C
tau = 150;
eps = eps1 * (1 + 2 * alpha * cos(g * X));
absE = sqrt(UB ./ eps / eps0); %场的模
btau = 1 / sqrt(-1/4 -1/2 * alpha * g * tau +1/4 * exp(2 * alpha * g * tau));
k1 = omega / c / sqrt(eps1);
x0tau = (g * X - 3 * pi / 2) / btau;
phaseE = k1 * exp(-g * alpha * sin(g * X) * tau) .* sqrt(pi) / 2 / g * btau .* (1 + erf(x0tau));
E = absE .* exp(1i * phaseE);

%图D
tau = 150;
UD = zeros(5, length(X));
deltalist = [-0.05, -0.025, 0, 0.025, 0.05];

for i = 1:length(deltalist)
    delta = deltalist(i);
    cg = (1 + delta) * c / sqrt(eps1);
    g = Omega / cg; %调制空间频率
    X = linspace(0, 2 * pi / g, 5000);
    UD(i, :) = exp(-g * (2 + delta) * alpha * sin(g * X) * tau - (delta + alpha * cos(g * X)) .^ 2/2 / alpha ^ 2 * (-1 - 2 * alpha * g * tau + exp(2 * alpha * g * tau)) ...
        +delta * (delta + alpha * cos(g * X)) * 2 / alpha ^ 2 * (-1 - alpha * g * tau + exp(alpha * g * tau)));

end

%% 绘图
figure(1)
plot(g * X / 2 / pi, UA, 'r-', 'LineWidth', 2)
xlim([0 1])
xlabel('gX/2\pi')
ylabel('U(X)')
grid on
title('能量密度 \tau=37.5')

figure(2)
plot(g * X / 2 / pi, UB, 'r-', 'LineWidth', 2)
xlim([0 1])
xlabel('gX/2\pi')
ylabel('U(X)')
grid on
title('能量密度 \tau=150')

figure(3)
plot(g * X / 2 / pi, real(E), 'r-', 'LineWidth', 2)
xlim([0.65, 0.85])
xlabel('gX/2\pi')
ylabel('Re(E(X))')
grid on
title('E(X)的幅度')

figure(4)
plot(g * X / 2 / pi, UD, 'LineWidth', 2)
xlim([0.5 1])
xlabel('gX/2\pi')
ylabel('U(X)')
grid on
legend('\delta=-0.05', '\delta=-0.025', '\delta=0', '\delta=0.025', '\delta=0.05')
title('不同偏差下的能量密度')
