%% 计算...的色散曲面 使用的方法是传输矩阵法

%初始化
clear
clc
close all

%电路参数
LR = 2.5 * 1e-9; %右手a长度电感
CR = 1e-12; %右手a长度电容
CV = 1e-12;

%扫频范围
fMin = 0;
fMax = 1e10;
f = fMin:fMax / 10000:fMax;
Omega = 2 * pi * fMax;
omega = 2 * pi * f;

%增量模型
Z = 1i * omega * LR;
Y = 1i * omega * CR;
C = 1 ./ (1i * omega * CV);
Zc = sqrt(LR / CR);
omega0 = sqrt(1 / (LR * CR));
N = 3; %y方向的层数，此处有三条线，所以是3
M = 12; %x方向的层数

%传输矩阵
T = zeros(2 * N, 2 * N, length(omega));
Th = zeros(2 * N, 2 * N, length(omega));
Tv = zeros(2 * N, 2 * N, length(omega));
Ttol = zeros(2 * N, 2 * N, length(omega));
X = zeros(N, N, length(omega));

%填充传输矩阵
for i = 1:length(omega)
    Th(:, :, i) = [eye(N), Z(i) / 2 * eye(N);
                   zeros(N), eye(N)];
    Ya = Y(i) + 1 / C(i);
    Yb = -1 / C(i);
    Yc = Y(i) + 2 / C(i);

    for k = 1:N

        for j = 1:N

            if k == j
                %主对角线
                if k == 1 || k == N
                    X(k, j, i) = Ya;
                else
                    X(k, j, i) = Yc;
                end

            end

            if abs(k - j) == 1
                X(k, j, i) = Yb;
            end

        end

    end

    Tv(:, :, i) = [
                   eye(N), zeros(N);
                   X(:, :, i), eye(N)
                   ];
    T(:, :, i) = Th(:, :, i) * Tv(:, :, i) * Th(:, :, i);
    Ttol(:, :, i) = T(:, :, i) ^ M;
end

Atol = Ttol(1:N, 1:N, :);
Btol = Ttol(1:N, N + 1:2 * N, :);
Ctol = Ttol(N + 1:2 * N, 1:N, :);
Dtol = Ttol(N + 1:2 * N, N + 1:2 * N, :);

%Z矩阵
Ztol = zeros(2 * N, 2 * N, length(omega));
Stol = zeros(2 * N, 2 * N, length(omega));

for i = 1:length(omega)
    Ztol(:, :, i) = [Atol(:, :, i) * Ctol(:, :, i) ^ -1, Ctol(:, :, i) ^ -1;
                     Ctol(:, :, i) ^ -1, Dtol(:, :, i) * Ctol(:, :, i) ^ -1];
    Stol(:, :, i) = (Ztol(:, :, i) / Zc + eye(2 * N)) * (Ztol(:, :, i) / Zc - eye(2 * N)) ^ -1;
end

S11 = squeeze(Stol(1, 1, :));
S41 = squeeze(Stol(4, 1, :));
S52 = squeeze(Stol(5, 2, :));
S63 = squeeze(Stol(6, 3, :));
S41_phi = angle(S41) .* 180 / pi;
S52_phi = angle(S52) .* 180 / pi;
S63_phi = angle(S63) .* 180 / pi;
S11_dB = 20 * log10(abs(S11));
S41_dB = 20 * log10(abs(S41));
S52_dB = 20 * log10(abs(S52));
S63_dB = 20 * log10(abs(S63));

%% 画图
h1 = figure;

subplot(4, 2, 1:2)
plot(f, S11_dB, '.', 'Color', 'b', 'MarkerSize', 1)
xlim([0 fMax])
ylim([-40 10])
yticks(-40:10:10)
title('|S11|')
grid on
subplot(4, 2, 3)
plot(f, S41_dB, '.', 'Color', 'b', 'MarkerSize', 1)
xlim([0 fMax])
ylim([-40 10])
yticks(-40:10:10)
title('|S41|')
grid on
subplot(4, 2, 4)
plot(f, S41_phi, '-', 'Color', 'b')
xlim([0 fMax])
ylim([-200 200])
title('S41 phase')
grid on
subplot(4, 2, 5)
plot(f, S52_dB, '.', 'Color', 'b', 'MarkerSize', 1)
xlim([0 fMax])
ylim([-40 10])
yticks(-40:10:10)
title('|S52|')
grid on
subplot(4, 2, 6)
plot(f, S52_phi, '-', 'Color', 'b')
xlim([0 fMax])
ylim([-200 200])
title('S52 phase')
grid on
subplot(4, 2, 7)
plot(f, S63_dB, '.', 'Color', 'b', 'MarkerSize', 1)
xlim([0 fMax])
ylim([-40 10])
yticks(-40:10:10)
title('|S63|')
grid on
subplot(4, 2, 8)
plot(f, S63_phi, '-', 'Color', 'b')
xlim([0 fMax])
ylim([-200 200])
title('S63 phase')
grid on
