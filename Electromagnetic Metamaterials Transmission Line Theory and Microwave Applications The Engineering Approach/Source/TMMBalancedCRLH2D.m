%% 验证TMM的可行性

%初始化
clear
clc
close all

%电路参数
LR = 2.5 * 1e-9; %右手单位长度电感
CR = 1e-12; %右手单位长度电感
LL = 2.5 * 1e-9; %右手单位长度电感
CL = 1e-12; %右手单位长度电感
Zth = 50;

%扫频范围
fMax = 1e10;
f = 0:fMax / 1000:fMax;
Omega = 2 * pi * fMax;
omega = 0:Omega / 1000:Omega;

Z = 1i * omega * LR + 1 ./ (1i * omega * CL);
Y = 1i * omega * CR + 1 ./ (1i * omega * LL);

Zc = sqrt(Z ./ Y);
omega0 = (LR * CR * LL * CL) ^ (-1/4);

N = 12; %层数

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
    Ya = Y(i) + 1 / Z(i);
    Yb = -1 / Z(i);
    Yc = Y(i) + 2 / Z(i);

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
    Ttol(:, :, i) = T(:, :, i) ^ N;
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
    Stol(:, :, i) = (Ztol(:, :, i) / Zc(i) + eye(2 * N)) * (Ztol(:, :, i) / Zc(i) - eye(2 * N)) ^ -1;
end

%% fig4.11
h1 = figure;
S11 = squeeze(Stol(1, 1, :));
S11_dB = 20 * log10(abs(S11));
S39 = squeeze(Stol(3, 9, :));
S39_dB = 20 * log10(abs(S39));
S66 = squeeze(Stol(6, 6, :));
S66_dB = 20 * log10(abs(S66));
S624 = squeeze(Stol(6, 24, :));
S624_dB = 20 * log10(abs(S624));

subplot(2, 2, 1)
plot(f, S11_dB, '-', 'Color', 'b')
ylim([-40 10])
xlabel('f')
ylabel('dB')
title('|S11|')
subplot(2, 2, 2)
plot(f, S39_dB, '-', 'Color', 'b')
ylim([-40 10])
xlabel('f')
ylabel('dB')
title('|S39|')
subplot(2, 2, 3)
plot(f, S66_dB, '-', 'Color', 'b')
ylim([-40 10])
xlabel('f')
ylabel('dB')
title('|S66|')
subplot(2, 2, 4)
plot(f, S624_dB, '-', 'Color', 'b')
ylim([-40 10])
xlabel('f')
ylabel('dB')
title('|S624|')
