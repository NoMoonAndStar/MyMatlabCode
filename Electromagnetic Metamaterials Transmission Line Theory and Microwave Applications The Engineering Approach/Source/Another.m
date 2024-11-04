%% 计算...的色散曲面

%初始化
clear
clc
close all

%电路参数
LR = 2.5 * 1e-9; %右手a/2长度电感
CR = 1e-12; %右手a长度电容
CV = 1e-12;

%扫频范围
Omega = 5e10;
omega = 0:Omega / 1000:Omega;

%增量模型
Z = 1i * omega * LR;
Y = 1i * omega * CR;
C = 1 ./ (1i * omega * CV);

%% 色散曲线
%gamma->X
kpx = acos(1 + Y .* Z / 2);
kpy = 0;
%这里把虚部为0的点去掉
idx = imag(kpx) == 0;
kpx = kpx(idx);
omega1 = omega(idx);

%因为Beta是向量，所以取范数
Beta1 = sqrt(real(kpy) .^ 2 + real(kpx) .^ 2);

h1 = figure;
subplot(1, 3, 1)
plot(Beta1, omega1, '.', 'Color', 'b')
ylim([0 Omega])
xlabel('Beta*p')
ylabel('omega')
title('Dispersion relation from gamma to X')

%X->M
kpx = pi;
kpy = acos(1 + C .* Y + 4 * C ./ Z);
idx = imag(kpy) == 0;
kpy = kpy(idx);
omega2 = omega(idx);
Beta2 = sqrt(real(kpy) .^ 2 + real(kpx) .^ 2);

subplot(1, 3, 2)
plot(Beta2, omega2, '.', 'Color', 'b')
ylim([0 Omega])
xlabel('Beta*p')
ylabel('omega')
title('Dispersion relation from X to M')

%M->gamma
kpx = acos(1 + C .* Y .* Z ./ (2 * C + Z));
kpy = acos(1 + C .* Y .* Z ./ (2 * C + Z));
idx = imag(kpx) == 0;
kpx = kpx(idx);
idx = imag(kpy) == 0;
kpy = kpy(idx);
omega3 = omega(idx);
Beta3 = sqrt(real(kpy) .^ 2 + real(kpx) .^ 2);

subplot(1, 3, 3)
plot(Beta3, omega3, '.', 'Color', 'b')
set(gca, 'XDir', 'reverse') %对X方向反转
ylim([0 Omega])
xlabel('Beta*p')
ylabel('omega')
title('Dispersion relation from M to gamma')

%% 色散曲面
m = 200;

KPX = nan(length(omega) - 1, m);
KPY = nan(length(omega) - 1, m);

for i = 2:length(omega) - 1
    [x_s, y_s] = solu2beyondequ(-2 * C(i), -Z(i), 2 * C(i) + Z(i) + Z(i) * C(i) * Y(i), m);

    if ~isempty(x_s)
        KPX(i, 1:length(x_s)) = x_s;
        KPY(i, 1:length(y_s)) = y_s;
    end

end

h4 = figure;
hold on

for i = 2:length(omega) - 1

    plot3(KPX(i, :), KPY(i, :), omega(i) * ones(1, m), 'b.')
end

xlabel('kx*px')
ylabel('ky*py')
zlabel('oemga')
title('Dispersion relation')
grid on
