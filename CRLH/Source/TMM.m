% 计算...的色散曲面 使用的方法是传输矩阵法

%初始化
clear
clc
close all

%电路参数
LR = 2.5 * 1e-9;
CR = 1e-12;
CV = 0.5e-12;

%主路扫频范围
f1Min = 0;
f1Max = 2e9;
f1 = f1Min:f1Max / 10000:f1Max;
Omega1 = 2 * pi * f1Max;
omega1 = 2 * pi * f1;

%调制路扫频范围
f2 = 300 * 1e6;
omega2 = 2 * pi * f2;

%主路阻抗
Z1 = 1i * omega1 * LR;
Y1 = 1i * omega1 * CR;

%调制路阻抗和电容阻抗
Z2 = 1i * omega2 * LR;
Y2 = 1i * omega2 * CR;
YC1 = 1i * omega1 * CV;
YC2 = 1i * omega2 * CV;
Zc = sqrt(LR / CR);

N = 3; %y方向的层数，此处有三条线，所以是3
M = 25; %x方向的层数
p = 10 * 1e-3;
L = p * M;

%传输矩阵
T = zeros(2 * N, 2 * N, length(omega1));
Th = zeros(2 * N, 2 * N, length(omega1));
Tv = zeros(2 * N, 2 * N, length(omega1));
Ttol = zeros(2 * N, 2 * N, length(omega1));
X = zeros(N, N, length(omega1));

%填充传输矩阵
for i = 1:length(omega1)
    Th(:, :, i) = [eye(N), diag([Z2 / 2, Z1(i) / 2, Z2 / 2]);
                                 zeros(N), eye(N)];
    X = [
         Y2 + YC2, -YC2, 0;
         -YC1(i), Y1(i) + 2 * YC1(i), -YC1(i);
         0, -YC2, Y2 + YC2
         ];

    Tv(:, :, i) = [
                   eye(N), zeros(N);
                   X, eye(N)
                   ];
    T(:, :, i) = Th(:, :, i) * Tv(:, :, i) * Th(:, :, i);
    Ttol(:, :, i) = T(:, :, i) ^ M;
end

Atol = Ttol(1:N, 1:N, :);
Btol = Ttol(1:N, N + 1:2 * N, :);
Ctol = Ttol(N + 1:2 * N, 1:N, :);
Dtol = Ttol(N + 1:2 * N, N + 1:2 * N, :);

%Z矩阵
Ztol = zeros(2 * N, 2 * N, length(omega1));
Stol = zeros(2 * N, 2 * N, length(omega1));

for i = 1:length(omega1)
    Ztol(:, :, i) = [Atol(:, :, i) * Ctol(:, :, i) ^ -1, Ctol(:, :, i) ^ -1;
                     Ctol(:, :, i) ^ -1, Dtol(:, :, i) * Ctol(:, :, i) ^ -1];
    Stol(:, :, i) = (Ztol(:, :, i) / Zc - eye(2 * N)) * (Ztol(:, :, i) / Zc + eye(2 * N)) ^ -1;
end

S11 = squeeze(Stol(1, 1, :));
S41 = squeeze(Stol(4, 1, :));
S52 = squeeze(Stol(5, 2, :));
S63 = squeeze(Stol(6, 3, :));
S41_phi = angle(S41) .* 180 / pi;
S52_phi = angle(S52) .* 180 / pi;
S63_phi = angle(S63) .* 180 / pi;
unwrapped_S41_phi = unwrap(S41_phi, -360);
unwrapped_S52_phi = unwrap(S52_phi, -360);
unwrapped_S63_phi = unwrap(S63_phi, -360);
beta41 = -unwrapped_S41_phi * pi / 180 / L;
beta52 = -unwrapped_S52_phi * pi / 180 / L;
beta63 = -unwrapped_S63_phi * pi / 180 / L;
theta41 = beta41 * p;
theta52 = beta52 * p;
theta63 = beta63 * p;
S11_dB = 20 * log10(abs(S11));
S41_dB = 20 * log10(abs(S41));
S52_dB = 20 * log10(abs(S52));
S63_dB = 20 * log10(abs(S63));

%% 画图
h1 = figure;

subplot(3, 4, 1)
plot(f1, S41_dB, 'b.')
xlim([0 f1Max])
ylim([-40 10])
yticks(-40:10:10)
title('|S41|')
grid on
subplot(3, 4, 2)
plot(f1, S41_phi, '-', 'Color', 'b')
xlim([0 f1Max])
ylim([-200 200])
title('S41 phase')
grid on
subplot(3, 4, 3)
plot(f1, unwrapped_S41_phi, 'r-')
grid on
subplot(3, 4, 4)
plot(theta41, f1, 'r-')
ylabel('f(Hz)')
xlabel('beta41*p')
grid on

subplot(3, 4, 5)
plot(f1, S52_dB, 'b.')
xlim([0 f1Max])
ylim([-40 10])
yticks(-40:10:10)
title('|S52|')
grid on
subplot(3, 4, 6)
plot(f1, S52_phi, '-', 'Color', 'b')
xlim([0 f1Max])
ylim([-200 200])
title('S52 phase')
grid on
subplot(3, 4, 7)
plot(f1, unwrapped_S52_phi, 'r-')
grid on
subplot(3, 4, 8)
plot(f1, omega1./beta52', 'r-')
ylabel('vp')
xlabel('f')
grid on

subplot(3, 4, 9)
plot(f1, S63_dB, 'b.')
xlim([0 f1Max])
ylim([-40 10])
yticks(-40:10:10)
title('|S63|')
grid on
subplot(3, 4, 10)
plot(f1, S63_phi, '-', 'Color', 'b')
xlim([0 f1Max])
ylim([-200 200])
title('S63 phase')
grid on
subplot(3, 4, 11)
plot(f1, unwrapped_S63_phi, 'r-')
grid on
subplot(3, 4, 12)
plot(theta63, f1, 'r-')
ylabel('f(Hz)')
xlabel('beta63*p')
grid on
