%% 初始化
clear
clc
close all

%% 参数定义
Z=50;
c=3e8;
a=0.2;
C=8.47*1e-12;
eps=1e-3;
m=500;
fmin=10e6; fmax=1e9; fnum=500;
f=linspace(fmin,fmax,fnum);
Omega=2*pi*f;

%% 解方程
%gamma->x
k2a=0; y=exp(-1i*k2a);
Result1=NaN(m*fnum,3);
count=1;
for omega = Omega
    k=omega/c;
    Ya=1i*omega*C;
    Yb=-1i*omega*C;
    Yc=2i*omega*C;
    X=[Ya Yb 0;Yb Yc Yb;0 Yb Ya];
    Th=[cos(k*a/2)*eye(3) 1i*Z*sin(k*a/2)*eye(3);1i*sin(k*a/2)*eye(3)/Z cos(k*a/2)*eye(3)];
    Tv=[eye(3) zeros(3);X eye(3)];
    T=Th*Tv*Th;
    for k1a = linspace(0,pi,m)
        x1=exp(-1i*k1a); x2=exp(1i*k1a);
        E=diag([x1,y,x1,x1,y,x1]);
        if abs(det(T-E)) < eps
            Result1(count,:)=[k1a,omega,abs(det(T-E))];
            count=count+1;
        end
    end
end

%% 绘图
plot(Result1(:,1),Result1(:,2)/2/pi,'b.')
xlim([0 pi])
ylim([fmin,fmax])
xlabel('k1a')
ylabel('k2a')
zlabel('Frequency')