%% 计算二维的CRLH的色散曲面

%初始化
clear
clc
close all

%四个电路参数
LR = 2.5 * 1e-9; %右手单位长度电感
CR = 1e-12; %右手单位长度电感
LL = 2.5 * 1e-9; %右手单位长度电感
CL = 1e-12; %右手单位长度电感

%扫频范围
Omega = 1e11;
omega = 0:Omega / 20:Omega;

%晶格间距
px = 10 * 1e-3;
py = 20 * 1e-3;

%电路
Z = 1i * omega * LR / 2 + 1 ./ (1i * 2 * omega * CL);
Y = 1i * omega * CR + 1 ./ (1i * omega * LL);
x = -Z .* Y;

%% fig4.5(a)
%gamma->X
kpx = acos(1 - x / 2);
kpy = 0;
%这里把虚部为0的点去掉
idx = imag(kpx) == 0;
kpx = kpx(idx);
omega1 = omega(idx);

%因为Beta是向量，所以取范数
Beta = sqrt(real(kpy) .^ 2 + real(kpx) .^ 2);

h1 = figure;
subplot(1, 3, 1)
plot(Beta, omega1, '.', 'Color', 'b')
ylim([0 Omega])
xlabel('Beta*p')
ylabel('omega')
title('Dispersion relation from gamma to X')

%X->M
kpx = pi;
kpy = acos(3 - x / 2);
idx = imag(kpy) == 0;
kpy = kpy(idx);
omega2 = omega(idx);

Beta = sqrt(real(kpy) .^ 2 + real(kpx) .^ 2);

subplot(1, 3, 2)
plot(Beta, omega2, '.', 'Color', 'b')
ylim([0 Omega])
xlabel('Beta*p')
ylabel('omega')
title('Dispersion relation from X to M')

%M->gamma
kpx = acos(1 - x / 4);
kpy = acos(1 - x / 4);
idx = imag(kpx) == 0;
kpx = kpx(idx);
idx = imag(kpy) == 0;
kpy = kpy(idx);

omega3 = omega(idx);

Beta = sqrt(real(kpy) .^ 2 + real(kpx) .^ 2);
%这里反转Beta是为了便于画图，没有实际意义
Beta = -Beta + max(Beta);

subplot(1, 3, 3)
plot(Beta, omega3, '.', 'Color', 'b')
ylim([0 Omega])
xlabel('Beta*p')
ylabel('omega')
title('Dispersion relation from M to gamma')

%% fig4.5(b)
kpx = linspace(-pi, pi, 100);
kpy = linspace(-pi, pi, 100);

% 创建网格
[KPX, KPY] = meshgrid(kpx, kpy);
Z = cos(KPX) + cos(KPY);
k = 2 - x / 2;
%取实数解
idx = k >= -2 & k <= 2;
k = k(idx);
omega4 = omega(idx);
eps = 1e-3;
results_kpx = cell(1, length(k));
results_kpy = cell(1, length(k));

for i = 1:length(k)
    xy_pairs = find(abs(Z - k(i)) < eps);
    [x_sol, y_sol] = ind2sub(size(Z), xy_pairs);
    results_kpx{i} = kpx(x_sol);
    results_kpy{i} = kpy(y_sol);
end

max_length = max(cellfun(@numel, results_kpx));

% 创建一个大数组并用NaN填充
big_array_kpx = NaN(length(k), max_length);
big_array_kpy = NaN(length(k), max_length);

% 填充数组
for i = 1:length(k)
    len = numel(results_kpx{i});
    big_array_kpx(i, 1:len) = results_kpx{i};
    big_array_kpy(i, 1:len) = results_kpy{i};
    %现在big_array每一行对应一个频率的kpx和kpy的值
end

% 将数组转换为列向量
X = big_array_kpx';
Y = big_array_kpy';
X = X(:)';
Y = Y(:)';
% 将频率转换为对角矩阵
Z = Remap(omega4, max_length);
Z(Z == 0) = NaN;

h2 = figure;
plot3(X, Y, Z, 'b.')
xlabel('kx*px')
ylabel('ky*py')
zlabel('omega')
title('Dispersion relation of 2D CRLH')
grid on
