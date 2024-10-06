clear
clc
close all
tic
%% 导入CST中的电场和磁场（时谐场）
% 加载电场数据
Ed = importdata("D:\CST\S-TheroyTEM\Export\3d\e-field (f=1) [pw].csv", ';', 2).data;
% 加载磁场数据
Hd = importdata("D:\CST\S-TheroyTEM\Export\3d\h-field (f=1) [pw].csv", ';', 2).data;
% 获取数据长度
N = length(Ed(:, 1));

% 初始化结构体存储点的信息
points = struct('x', [], 'y', [], 'z', [], 'Ex', [], 'Ey', [], 'Ez', [], 'Hx', [], 'Hy', [], 'Hz', [], 'Sx', [], 'Sy', [], 'Sz', []);
points = repmat(points, N, 1);

% 遍历所有数据点
for i = 1:N
    points(i).x = Ed(i, 1) * 1e-3; % x坐标
    points(i).y = Ed(i, 2) * 1e-3; % y坐标
    points(i).z = Ed(i, 3) * 1e-3; % z坐标
    points(i).Ex = complex(Ed(i, 4), Ed(i, 5)); % 复数形式的Ex
    points(i).Ey = complex(Ed(i, 6), Ed(i, 7)); % 复数形式的Ey
    points(i).Ez = complex(Ed(i, 8), Ed(i, 9)); % 复数形式的Ez
    points(i).Hx = complex(Hd(i, 4), Hd(i, 5)); % 复数形式的Hx
    points(i).Hy = complex(Hd(i, 6), Hd(i, 7)); % 复数形式的Hy
    points(i).Hz = complex(Hd(i, 8), Hd(i, 9)); % 复数形式的Hz
    points(i).Sx = 1/2 * (points(i).Ey * conj(points(i).Hz) - points(i).Ez * conj(points(i).Hy)); % 复坡印廷矢量
    points(i).Sy = 1/2 * (points(i).Ez * conj(points(i).Hx) - points(i).Ex * conj(points(i).Hz));
    points(i).Sz = 1/2 * (points(i).Ex * conj(points(i).Hy) - points(i).Ey * conj(points(i).Hx));
end

% 设置介质参数
epsr = 5; % 相对介电常数=' .
mur = 1; % 相对磁导率
eps0 = 8.854e-12; % 真空介电常数
mu0 = 4 * pi * 1e-7; % 真空磁导率
tanD = 0; % 损耗角正切
epsi = tanD * eps0 * epsr;
mui = 0;
f = 1 * 1e9; % 频率
dx = 1 * 1e-3; % 积分区域尺寸
dy = 10 * 1e-3;
dz = 1 * 1e-3;

%% 求右边
% 储能
w_store = zeros(N, 1);

for i = 1:N
    w_store(i) = 1/2 * 2 * pi * f * ((mu0 * mur * sum([points(i).Hx, points(i).Hy, points(i).Hz] .* conj([points(i).Hx, points(i).Hy, points(i).Hz]))) - eps0 * epsr * sum([points(i).Ex, points(i).Ey, points(i).Ez] .* conj([points(i).Ex, points(i).Ey, points(i).Ez])));
end

iright = 1i * mean(w_store) * dx * dy * dz;

% 损耗
w_loss = zeros(N, 1);

for i = 1:N
    w_loss(i) = 1/2 * 2 * pi * f * ((mui * sum([points(i).Hx, points(i).Hy, points(i).Hz] .* conj([points(i).Hx, points(i).Hy, points(i).Hz]))) + epsi * sum([points(i).Ex, points(i).Ey, points(i).Ez] .* conj([points(i).Ex, points(i).Ey, points(i).Ez])));
end

rright = mean(w_loss) * dx * dy * dz;
right = iright + rright;

%% 求左边：坡印廷矢量的面积分

Sint_x = [];
Sint_y = [];
Sint_z = [];

for i = 1:N

    if points(i).y == 0
        Sint_y = [Sint_y, -points(i).Sy];
    elseif points(i).y == dy
        Sint_y = [Sint_y, points(i).Sy];
    elseif points(i).x == -dx / 2
        Sint_x = [Sint_x, -points(i).Sx];
    elseif points(i).x == dx / 2
        Sint_x = [Sint_x, points(i).Sx];
    elseif points(i).z == 0
        Sint_z = [Sint_z, -points(i).Sz];
    elseif points(i).z == dz
        Sint_z = [Sint_z, points(i).Sz];
    end

end

left =- 2 * (mean(Sint_x) * dy * dz + mean(Sint_y) * dz * dx + mean(Sint_z) * dx * dy);

%% 对比结果
table(left, right)
toc
