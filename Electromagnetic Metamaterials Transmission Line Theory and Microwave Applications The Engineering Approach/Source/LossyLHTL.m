%%有损耗的LH TL的传播常数和特性阻抗

%%初始化
clear
clc
close all

%电路参数
CL = 1e-6; %左手倍长度电容
LL = 1e-6; %左手倍长度电感
omegaL = sqrt(1 / (LL * CL));

%损耗
R = [0; 1; 10];
G = [0; 1e-6; 1e-5];

Omega = 3e7;
omega = 0:Omega / 20000:Omega; %扫频范围

%%fig3.11(a)
gammaL = -1i * sqrt(1 - R .* G * (omega / omegaL) .^ 2 + 1i * (CL * R + LL * G) * omega) ./ (omega / omegaL); %传播常数
alpha = real(gammaL); %衰减常数
beta = imag(gammaL); %相位常数

h1 = figure;
hold on
COLOR = ['r', 'b', 'g'];

for i = 1:3
    plot(beta(i, :), omega, '.', 'Color', COLOR(i))
end

for i = 1:3
    plot(alpha(i, 2:end), omega(2:end), '--', 'Color', COLOR(i))
end

plot([0 0], [0, Omega], "LineWidth", 0.5, "Color", "k")
xlim([-1 1])
xlabel('beta/alpha')
ylabel('omega')
title('Dispersion relation')
legend('loss-less', 'weakly lossy', 'strongly lossy')

%%fig3.11(b)
ZL = sqrt(LL / CL);
Zc = sqrt((R - 1i ./ (omega * CL)) ./ (G - 1i ./ (omega * LL)));

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

xlabel('Zc/ZL')
ylabel('omega')
title('Characteristic impedance')
legend('loss-less', 'weakly lossy', 'strongly lossy')
