%% 复现3.1和3.2的图
%初始化
clear
clc
close all

%四个电路参数
LR = 1e-8; %右手单位长度电感
CL = 1e-6; %左手倍长度电容
CR = 2e-8; %右手单位长度电感
LL = 5e-6; %左手倍长度电感

%中间量
omegaR = sqrt(1 / (LR * CR));
omegaL = sqrt(1 / (LL * CL));
k = LR * CL + LL * CR;
omegase = sqrt(1 / (LR * CL));
omegash = sqrt(1 / (LL * CR));

Omega = 3e7;
omega = 0:Omega / 20000:Omega; %扫频范围

%损耗
R = [0; 1; 10];
G = [0; 1e-6; 1e-5];

Z = R + 1i * (omega * LR - 1 ./ (omega * CL));
Y = G + 1i * (omega * CR - 1 ./ (omega * LL));

%%fig3.13(a)
gamma = sqrt(Z .* Y);
alpha = real(gamma);
beta = imag(gamma);

idx = omega < min(omegash, omegase);
beta(1, 1:sum(idx)) = -beta(1, 1:sum(idx));

h1 = figure;
hold on
COLOR = ['r', 'b', 'g'];

for i = 1:3
    plot(beta(i, :), omega, '.', 'Color', COLOR(i))
end

for i = 1:3
    plot(alpha(i, 2:end), omega(2:end), '--', 'Color', COLOR(i))
end

plot([0 0], [0, Omega], "LineWidth", 1, "Color", "k")
xlim([-1 1])
xlabel('beta/alpha')
ylabel('omega')
title('Dispersion relation of Lossy Unbalanced CRLH')
legend('loss-less', 'weakly lossy', 'strongly lossy')

%%fig3.13(b)
Zc = sqrt(Z ./ Y);
ZL = sqrt(LL / CL);
Rc = real(Zc);
Xc = imag(Zc);

h2 = figure;
hold on

for i = 1:3
    plot(Rc(i, :) / ZL, omega, '.', 'Color', COLOR(i))
end

for i = 1:3
    plot(Xc(i, :) / ZL, omega, '--', 'Color', COLOR(i))
end

plot([0 0], [0, Omega], "LineWidth", 1, "Color", "k")
xlim([-5 15])
xlabel('Zc/ZL')
ylabel('omega')
title('Characteristic impedance of Lossy Unbalanced CRLH')
legend('loss-less', 'weakly lossy', 'strongly lossy')