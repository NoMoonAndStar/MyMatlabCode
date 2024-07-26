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
Omega = 5e10;
omega = 0:Omega / 10000:Omega;

%增量模型
Z = 1i * omega * LR;
Y = 1i * omega * CR;
C = 1 ./ (1i * omega * CV);
Zc = sqrt(LR / CR);
omega0 = sqrt(1 / (LR * CR));
N = 25; %层数

%传输矩阵
T = zeros(6, 6, length(omega));
Th = zeros(6, 6, length(omega));
Tv = zeros(6, 6, length(omega));
Ttol = zeros(6, 6, length(omega));
X = zeros(3, 3, length(omega));

%填充传输矩阵
for i = 1:length(omega)
    Th(:, :, i) = [diag(ones(1, 3)), Z(i) / 2 * diag(ones(1, 3));
                   zeros(3, 3), diag(ones(1, 3))];
    Ya = Y(i) + 1 / C(i);
    Yb = -1 / C(i);
    Yc = Y(i) + 2 / C(i);
    X(:, :, i) = [
                  Ya, Yb, 0;
                  Yb, Yc, Yb;
                  0, Yb, Ya
                  ];
    Tv(:, :, i) = [
                   diag(ones(1, 3)), zeros(3, 3);
                   X(:, :, i), diag(ones(1, 3))
                   ];
    T(:, :, i) = Th(:, :, i) * Tv(:, :, i) * Th(:, :, i);
    Ttol(:, :, i) = T(:, :, i) ^ N;
end

Atol = Ttol(1:3, 1:3, :);
Btol = Ttol(4:6, 1:3, :);
Ctol = Ttol(1:3, 4:6, :);
Dtol = Ttol(4:6, 4:6, :);

%Z矩阵
Ztol = zeros(6, 6, length(omega));
Stol = zeros(6, 6, length(omega));

for i = 1:length(omega)
    Ztol(:, :, i) = [Atol(:, :, i) * Ctol(:, :, i) ^ -1, Ctol(:, :, i) ^ -1;
                     Ctol(:, :, i) ^ -1, Ctol(:, :, i) ^ -1 * Dtol(:, :, i)];
    Stol(:, :, i) = (Ztol(:, :, i) / Zc + diag(ones(1, 6))) * (Ztol(:, :, i) / Zc - diag(ones(1, 6))) ^ -1;
end

S11 = squeeze(Stol(1, 1, :));
S41 = squeeze(Stol(4, 1, :));
S52 = squeeze(Stol(5, 2, :));
S41_phi = angle(S41) .* 180 / pi;
S52_phi = angle(S52) .* 180 / pi;
idx_omega0 = find(omega == omega0) + 1;
unwrapped_S41_phi = unwrap(S41_phi, -360);
unwrapped_S41_phi = unwrapped_S41_phi - unwrapped_S41_phi(idx_omega0);
unwrapped_S52_phi = unwrap(S52_phi, -360);
unwrapped_S52_phi = unwrapped_S52_phi - unwrapped_S52_phi(idx_omega0);
S11_dB = 20 * log10(abs(S11));
S41_dB = 20 * log10(abs(S41));
S52_dB = 20 * log10(abs(S52));

h1 = figure;
subplot(3, 2, 1:2)
plot(omega, S11_dB, '-', 'Color', 'b')
title('|S11|')
subplot(3, 2, 3)
plot(omega, S41_dB, '-', 'Color', 'b')
title('|S41|')
subplot(3, 2, 4)
plot(omega, unwrapped_S41_phi, '-', 'Color', 'b')
title('S41 phase')
subplot(3, 2, 5)
plot(omega, S52_dB, '-', 'Color', 'b')
title('|S52|')
subplot(3, 2, 6)
plot(omega, unwrapped_S52_phi, '-', 'Color', 'b')
title('S52 phase')
