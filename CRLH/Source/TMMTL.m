%% 计算...的色散曲面 使用的方法是传输矩阵法

%初始化
clear
clc
close all

%电路参数
E = 48 * pi / 180;
F = 1e9;
c = 3e8;
CV = 8.47e-12;
a = c * E / F / 2 / pi;

%扫频范围
fMin = 0;
fMax = 1e9;
f = fMin:fMax / 100000:fMax;
Omega = 2 * pi * fMax;
omega = 2 * pi * f;

%增量模型
Z0 = 50;
Y0 = 1 / Z0;
C = 1 ./ (1i * omega * CV);
N = 3; %y方向的层数，此处有三条线，所以是3
M = 64; %x方向的层数
L = a * M;

%传输矩阵
T = zeros(2 * N, 2 * N, length(omega));
Th = zeros(2 * N, 2 * N, length(omega));
Tv = zeros(2 * N, 2 * N, length(omega));
Ttol = zeros(2 * N, 2 * N, length(omega));
X = zeros(N, N, length(omega));

%填充传输矩阵
for i = 1:length(omega)
    % ?为什么这里用f而不是Omega
    theta = E / F * f(i);
    %这里传输矩阵是未归一化的
    Th(:, :, i) = [cos(theta / 2) * eye(3), 1i * Z0 * sin(theta / 2) * eye(3);
                   1i * Y0 * sin(theta / 2) * eye(3), cos(theta / 2) * eye(3)
                   ];
    Ya = 1 / C(i);
    Yb = -1 / C(i);
    Yc = 2 / C(i);

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
    Stol(:, :, i) = (Ztol(:, :, i) / Z0 - eye(2 * N)) * (Ztol(:, :, i) / Z0 + eye(2 * N)) ^ -1;
end

Z11 = squeeze(Ztol(1, 1, :));
Z22 = squeeze(Ztol(2, 2, :));
Z44 = squeeze(Ztol(1, 1, :));
Z55 = squeeze(Ztol(1, 1, :));
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
theta41 = beta41 * a;
theta52 = beta52 * a;
theta63 = beta63 * a;
S11_dB = 20 * log10(abs(S11));
S41_dB = 20 * log10(abs(S41));
S52_dB = 20 * log10(abs(S52));
S63_dB = 20 * log10(abs(S63));

%% 画图
h1 = figure;
plot(f, S41_dB, 'b-')
ylim([-40 10])
xlim([0 fMax])
yticks(-40:10:10)
title('S41（调制路）的幅度')
grid on

h2 = figure;
plot(f, S41_phi, '-', 'Color', 'b')
xlim([0 fMax])
ylim([-200 200])
title('S41（调制路）的相位（未展开）')
grid on

h3 = figure;
plot(f, unwrapped_S41_phi, 'r-')
title('S41（调制路）的相位（展开后）')
grid on

h4 = figure;
plot(beta41, f, 'r-')
ylabel('f(Hz)')
xlabel('β2')
title('调制路的色散关系')
grid on

h5 = figure;
plot(f, S52_dB, 'b-')
xlim([0 fMax])
ylim([-40 10])
yticks(-40:10:10)
title('S52（主路）的幅度')
grid on

h6 = figure;
plot(f, S52_phi, '-', 'Color', 'b')
xlim([0 fMax])
ylim([-200 200])
title('S52（主路）的相位（未展开）')
grid on

h7 = figure;
plot(f, unwrapped_S52_phi, 'r-')
title('S52（主路）的相位（展开后）')
grid on

h8 = figure;
plot(beta52, f, 'r-')
ylabel('f(Hz)')
xlabel('β1')
title('主路的色散关系')
grid on

% 相速度
vp1 = omega' ./ beta52;
vp2 = omega' ./ beta41;

h9 = figure;
hold on
plot(f, vp1, 'b-')
plot(f, vp2, 'r-')
grid on
legend("主路", "调制路")
