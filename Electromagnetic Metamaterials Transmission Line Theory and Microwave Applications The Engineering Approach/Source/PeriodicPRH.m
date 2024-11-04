%%周期 对称 无耗 PRH的特性

%初始化
clear
clc
close all

LR = 2.5 * 1e-9; %右手单位长度电感
CR = 1e-12; %右手单位长度电感

Omega = 1e11;
omega = 0:Omega / 10000:Omega; %扫频范围

p = 3.75e-4; %网格间距
N = 10; %单元个数
L = p * N; %网格长度

Z = 1i * omega * LR / 2;
Y = 1i * omega * CR;

Zc = sqrt(Z ./ Y);

A = zeros(2, 2, length(omega)); %单元传输矩阵
AN = zeros(2, 2, length(omega)); %总传输矩阵
SN = zeros(2, 2, length(omega)); %散射矩阵

for i = 1:length(omega)
    A(:, :, i) = [1 Z(i) / 2; 0 1] * [1 0; Y(i) 1] * [1 Z(i) / 2; 0 1];
    AN(:, :, i) = A(:, :, i) ^ N;
    SN(:, :, i) = 1 / (AN(1, 1, i) + AN(1, 2, i) / Zc(i) + AN(2, 1, i) * Zc(i) + AN(2, 2, i)) * ...
        [
     AN(1, 1, i) + AN(1, 2, i) / Zc(i) - AN(2, 1, i) * Zc(i) - AN(2, 2, i), 2 * (AN(1, 1, i) * AN(2, 2, i) - AN(1, 2, i) * AN(2, 1, i));
     2, -AN(1, 1, i) + AN(1, 2, i) / Zc(i) - AN(2, 1, i) * Zc(i) + AN(2, 2, i)
     ];
end

S11 = squeeze(SN(1, 1, :));
S21 = squeeze(SN(2, 1, :));

%%fig3.20(a)
h1 = figure;
hold on
S11_dB = 20 * log10(abs(S11));
S21_dB = 20 * log10(abs(S21));
plot(omega, S11_dB, '-', 'Color', 'b')
plot(omega, S21_dB, '-', 'Color', 'r')
ylim([-60 10])
xlabel('omega')
ylabel('magnitude')
title('S11/S21')
legend('S11', 'S21')

%%fig3.20(b)
h2 = figure;
hold on
S11_phi = angle(S11) .* 180 / pi;
S21_phi = angle(S21) .* 180 / pi;
plot(omega, S21_phi, '-', 'Color', 'b')
xlabel('omega')
ylabel('phase')
title('S21 phase')

%%fig3.20(c)
h3 = figure;
hold on
unwrappedS21_phi = unwrap(S21_phi, -360);
plot(omega, unwrappedS21_phi, '-', 'Color', 'b')
ylim([-2000, 0])
xlabel('omega')
ylabel('phase')
title('unwrapped S21 phase')

%%fig3.20(d)
beta = -unwrappedS21_phi / L; %相位常数
alpha = -log(abs(S21)) / L; %衰减常数
theta = beta * p; %电长度

h4 = figure;
plot(theta, omega, '.', 'color', 'b')
xlabel('theta')
ylabel('omega')
title('Dispersion relation')
