%% 用PRL的理论和JOSB的图对比
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
eps = eps0 * eps1;
mu0 = 4 * pi * 1e-7;
c = 3e8;
cg = c / sqrt(eps1);
g = Omega / cg; %调制空间频率

%% 计算
X = linspace(0, 2 * pi / g, 5000);
%图A
tau = 37.5;
t = tau / cg;
%根据Supplementary Material for ”Broadband Nonreciprocal Amplification in Luminal Metamaterials”中式(12)求能量密度
UA = exp(-2 * alpha * Omega * t * sin(g * X) - alpha ^ 2 * Omega ^ 2 * t ^ 2 * (cos(g * X)) .^ 2);

%图B
tau = 150;
t = tau / cg;
UB = exp(-2 * alpha * Omega * t * sin(g * X) - alpha ^ 2 * Omega ^ 2 * t ^ 2 * (cos(g * X)) .^ 2);

%图C
%根据Supplementary Material for ”Broadband Nonreciprocal Amplification in Luminal Metamaterials”中第二节求场
tau = 150;
t = tau / cg;
d = t * cg;
N = 100; %展开阶数
Nt = (-N:N)';
Nnum = 2 * N + 1;

%求真空里的特征值和特征向量
alpha = 0;
M = eigMatrix(Nt, g, omega, Omega, eps0, alpha);
[v0, k0] = eig(M);
k0 = diag(k0);
[k0, V] = sort(real(k0));
v0 = v0(:, V);

%求介质里的特征值和特征向量
alpha = 0.05;
M = eigMatrix(Nt, g, omega, Omega, eps0, alpha);
[v1, k1] = eig(M);
k1 = diag(k1);
[k1, V] = sort(real(k1));
v1 = v1(:, V);

%特征矩阵
Mvinc = v0(:, 1:Nnum);
Mvref = v0(:, Nnum + 1:2 * Nnum);
Mm = v1;

%辅助矩阵
evinc = [zeros(N, 1); 1; zeros(N, 1)];
P = diag(exp(1i * d * k1));
A = Mm * (Mm * P) ^ -1 * Mvinc;
B = -Mvref;
AB = [A, B];
C = [Mvinc, zeros(2 * Nnum, Nnum)] * [evinc; zeros(Nnum, 1)];

%解方程
eall = AB \ C;
etra = eall(1:Nnum);
EH = Mvinc * etra;

%合成场
k1 = omega / c / sqrt(eps1);

t = linspace(d - 2 * pi / g, d, 5000) / cg;
et = exp(1i * ((k1 + Nt * g) * d - (omega + Nt * Omega) * t));
ett = [et; et];
EHtra = EH .* ett;
Etras = sum(EHtra(1:Nnum, :));

E = abs(Etras) .^ 2;
E = Etras;

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
xlabel('gX/2\pi')
ylabel('Re(E(X))')
grid on
title('E(X)的幅度')
