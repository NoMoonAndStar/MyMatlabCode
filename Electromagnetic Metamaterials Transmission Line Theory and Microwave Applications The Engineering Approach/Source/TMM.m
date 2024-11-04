%% 计算...的色散曲面 使用的方法是传输矩阵法

%初始化
clear
clc
close all

%电路参数
LR = 2.5 * 1e-9; %右手a长度电感
CR = 1e-12; %右手a长度电容
CV = 0.5e-12;

%扫频范围
fMin = 0;
fMax = 2e9;
f = fMin:fMax / 1000:fMax;
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
p = 10 * 1e-3;
L = p * M;

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

subplot(2, 4, 1)
plot(f, S41_dB, 'b-')
xlim([0 fMax])
ylim([-40 10])
yticks(-40:10:10)
title('S41（调制路）的幅度')
grid on
subplot(2, 4, 2)
plot(f, S41_phi, '-', 'Color', 'b')
xlim([0 fMax])
ylim([-200 200])
title('S41（调制路）的相位（未展开）')
grid on
subplot(2, 4, 3)
plot(f, unwrapped_S41_phi, 'r-')
title('S41（调制路）的相位（展开后）')
grid on
subplot(2, 4, 4)
plot(theta41, f, 'r-')
ylabel('f(Hz)')
xlabel('β2p')
title('调制路的色散关系')
grid on

subplot(2, 4, 5)
plot(f, S52_dB, 'b-')
xlim([0 fMax])
ylim([-40 10])
yticks(-40:10:10)
title('S52（主路）的幅度')
grid on
subplot(2, 4, 6)
plot(f, S52_phi, '-', 'Color', 'b')
xlim([0 fMax])
ylim([-200 200])
title('S52（主路）的相位（未展开）')
grid on
subplot(2, 4, 7)
plot(f, unwrapped_S52_phi, 'r-')
title('S52（主路）的相位（展开后）')
grid on
subplot(2, 4, 8)
plot(theta52, f, 'r-')
ylabel('f(Hz)')
xlabel('β1p')
title('主路的色散关系')
grid on
