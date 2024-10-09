clear
clc
close all
tic
%% 设置介质参数
epsr = 5; % 相对介电常数=' .
mur = 1; % 相对磁导率
eps0 = 8.854e-12; % 真空介电常数
mu0 = 4 * pi * 1e-7; % 真空磁导率
tanD = 0.1; % 损耗角正切
epsi = tanD * epsr;
mui = 0;
f = 1 * 1e9; % 频率
dx = 1 * 1e-3; % 积分区域尺寸
dy = 10 * 1e-3;
dz = 1 * 1e-3;

%% 计算xm面上的面积分
% 加载电场数据
Ed = importdata("D:\CST\S-TheroyTEM\Export\3d\Exm_e-field (f=1) [pw].csv", ';', 2).data;
% 加载磁场数据
Hd = importdata("D:\CST\S-TheroyTEM\Export\3d\Hxm_h-field (f=1) [pw].csv", ';', 2).data;
% 获取数据长度
N = length(Ed(:, 1));

Ex = complex(Ed(:, 4), Ed(:, 5));
Ey = complex(Ed(:, 6), Ed(:, 7));
Ez = complex(Ed(:, 8), Ed(:, 9));
Hx = complex(Hd(:, 4), Hd(:, 5));
Hy = complex(Hd(:, 6), Hd(:, 7));
Hz = complex(Hd(:, 8), Hd(:, 9));
Sx = 1/2 * (Ey .* conj(Hz) - conj(Hy) .* Ez);
Sy = 1/2 * (Ez .* conj(Hx) - conj(Hz) .* Ex);
Sz = 1/2 * (Ex .* conj(Hy) - conj(Hx) .* Ex);

Sint_xm = -mean(Sx) * dy * dz;

%% 计算xm面上的面积分
% 加载电场数据
Ed = importdata("D:\CST\S-TheroyTEM\Export\3d\Exp_e-field (f=1) [pw].csv", ';', 2).data;
% 加载磁场数据
Hd = importdata("D:\CST\S-TheroyTEM\Export\3d\Hxp_h-field (f=1) [pw].csv", ';', 2).data;
% 获取数据长度
N = length(Ed(:, 1));

Ex = complex(Ed(:, 4), Ed(:, 5));
Ey = complex(Ed(:, 6), Ed(:, 7));
Ez = complex(Ed(:, 8), Ed(:, 9));
Hx = complex(Hd(:, 4), Hd(:, 5));
Hy = complex(Hd(:, 6), Hd(:, 7));
Hz = complex(Hd(:, 8), Hd(:, 9));
Sx = 1/2 * (Ey .* conj(Hz) - conj(Hy) .* Ez);
Sy = 1/2 * (Ez .* conj(Hx) - conj(Hz) .* Ex);
Sz = 1/2 * (Ex .* conj(Hy) - conj(Hx) .* Ex);

Sint_xp = mean(Sx) * dy * dz;

%% 计算ym面上的面积分
% 加载电场数据
Ed = importdata("D:\CST\S-TheroyTEM\Export\3d\Eym_e-field (f=1) [pw].csv", ';', 2).data;
% 加载磁场数据
Hd = importdata("D:\CST\S-TheroyTEM\Export\3d\Hym_h-field (f=1) [pw].csv", ';', 2).data;
% 获取数据长度
N = length(Ed(:, 1));

Ex = complex(Ed(:, 4), Ed(:, 5));
Ey = complex(Ed(:, 6), Ed(:, 7));
Ez = complex(Ed(:, 8), Ed(:, 9));
Hx = complex(Hd(:, 4), Hd(:, 5));
Hy = complex(Hd(:, 6), Hd(:, 7));
Hz = complex(Hd(:, 8), Hd(:, 9));
Sx = 1/2 * (Ey .* conj(Hz) - conj(Hy) .* Ez);
Sy = 1/2 * (Ez .* conj(Hx) - conj(Hz) .* Ex);
Sz = 1/2 * (Ex .* conj(Hy) - conj(Hx) .* Ex);

Sint_ym = -mean(Sy) * dx * dz;

%% 计算yp面上的面积分
% 加载电场数据
Ed = importdata("D:\CST\S-TheroyTEM\Export\3d\Eyp_e-field (f=1) [pw].csv", ';', 2).data;
% 加载磁场数据
Hd = importdata("D:\CST\S-TheroyTEM\Export\3d\Hyp_h-field (f=1) [pw].csv", ';', 2).data;
% 获取数据长度
N = length(Ed(:, 1));

Ex = complex(Ed(:, 4), Ed(:, 5));
Ey = complex(Ed(:, 6), Ed(:, 7));
Ez = complex(Ed(:, 8), Ed(:, 9));
Hx = complex(Hd(:, 4), Hd(:, 5));
Hy = complex(Hd(:, 6), Hd(:, 7));
Hz = complex(Hd(:, 8), Hd(:, 9));
Sx = 1/2 * (Ey .* conj(Hz) - conj(Hy) .* Ez);
Sy = 1/2 * (Ez .* conj(Hx) - conj(Hz) .* Ex);
Sz = 1/2 * (Ex .* conj(Hy) - conj(Hx) .* Ex);

Sint_yp = mean(Sy) * dx * dz;

%% 计算zm面上的面积分
% 加载电场数据
Ed = importdata("D:\CST\S-TheroyTEM\Export\3d\Ezm_e-field (f=1) [pw].csv", ';', 2).data;
% 加载磁场数据
Hd = importdata("D:\CST\S-TheroyTEM\Export\3d\Hzm_h-field (f=1) [pw].csv", ';', 2).data;
% 获取数据长度
N = length(Ed(:, 1));

Ex = complex(Ed(:, 4), Ed(:, 5));
Ey = complex(Ed(:, 6), Ed(:, 7));
Ez = complex(Ed(:, 8), Ed(:, 9));
Hx = complex(Hd(:, 4), Hd(:, 5));
Hy = complex(Hd(:, 6), Hd(:, 7));
Hz = complex(Hd(:, 8), Hd(:, 9));
Sx = 1/2 * (Ey .* conj(Hz) - conj(Hy) .* Ez);
Sy = 1/2 * (Ez .* conj(Hx) - conj(Hz) .* Ex);
Sz = 1/2 * (Ex .* conj(Hy) - conj(Hx) .* Ex);

Sint_zm = -mean(Sz) * dx * dy;

%% 计算zp面上的面积分
% 加载电场数据
Ed = importdata("D:\CST\S-TheroyTEM\Export\3d\Ezp_e-field (f=1) [pw].csv", ';', 2).data;
% 加载磁场数据
Hd = importdata("D:\CST\S-TheroyTEM\Export\3d\Hzp_h-field (f=1) [pw].csv", ';', 2).data;
% 获取数据长度
N = length(Ed(:, 1));

Ex = complex(Ed(:, 4), Ed(:, 5));
Ey = complex(Ed(:, 6), Ed(:, 7));
Ez = complex(Ed(:, 8), Ed(:, 9));
Hx = complex(Hd(:, 4), Hd(:, 5));
Hy = complex(Hd(:, 6), Hd(:, 7));
Hz = complex(Hd(:, 8), Hd(:, 9));
Sx = 1/2 * (Ey .* conj(Hz) - conj(Hy) .* Ez);
Sy = 1/2 * (Ez .* conj(Hx) - conj(Hz) .* Ex);
Sz = 1/2 * (Ex .* conj(Hy) - conj(Hx) .* Ex);

Sint_zp = mean(Sz) * dx * dy;

%% 计算长方体内的体积分
% 加载电场数据
Ed = importdata("D:\CST\S-TheroyTEM\Export\3d\Ev_e-field (f=1) [pw].csv", ';', 2).data;
% 加载磁场数据
Hd = importdata("D:\CST\S-TheroyTEM\Export\3d\Hv_h-field (f=1) [pw].csv", ';', 2).data;
% 获取数据长度
N = length(Ed(:, 1));

Ex = complex(Ed(:, 4), Ed(:, 5));
Ey = complex(Ed(:, 6), Ed(:, 7));
Ez = complex(Ed(:, 8), Ed(:, 9));
Hx = complex(Hd(:, 4), Hd(:, 5));
Hy = complex(Hd(:, 6), Hd(:, 7));
Hz = complex(Hd(:, 8), Hd(:, 9));
wes = 1/4 * epsr * eps0 * sum([Ex, Ey, Ez] .* conj([Ex, Ey, Ez]), 2);
wms = 1/4 * mur * mu0 * sum([Hx, Hy, Hz] .* conj([Hx, Hy, Hz]), 2);
wel = 1/2 * epsi * eps0 * sum([Ex, Ey, Ez] .* conj([Ex, Ey, Ez]), 2);
wml = 1/2 * mui * mu0 * sum([Hx, Hy, Hz] .* conj([Hx, Hy, Hz]), 2);

We = 2 * 2 * pi * f * mean(wes) * dx * dy * dz;
Wh = 2 * 2 * pi * f * mean(wms) * dx * dy * dz;
Wstore = 2 * 2 * pi * f * mean(wms - wes) * dx * dy * dz;
Wloss = 2 * pi * f * mean(wel + wml) * dx * dy * dz;

%% 整合结果
left =- (Sint_xm + Sint_xp + Sint_ym + Sint_yp + Sint_zm + Sint_zp);
right = complex(Wloss, Wstore);

table(We,Wh)
table(left, right)

