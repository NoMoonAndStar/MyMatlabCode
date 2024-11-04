% 数值解久期方程
clear
clc

%% 预定义矩阵
tic
fMin=0; fMax=1e9;
KX=[]; KY=[]; F=[];
c=3e8;
Z01=96.23;
Z02=27.7249;
l1=0.0065;
l2=0.0065;
C=8.47*1e-12;
m=100;
n=100;
eps=0.1;
kx=linspace(-pi,pi,m);
ky=linspace(-pi,pi,m);
omega=linspace(fMin,2*pi*fMax,n);

%% 计算
parpool('Processes', 6) %进程池
corenum=6;
if isempty(gcp('nocreate'))
    parpool(corenum)
else
    delete(gcp('nocreate'))
    parpool(corenum)
end

parfor o = 1:n
    for i=1:m
        for j=1:m
            x1=exp(-1i*kx(i)); x2=exp(1i*kx(i)); y=exp(-1i*ky(j));
            k=omega(o)/c;
            Ya=1i*omega(o)*C;
            Yb=-1i*omega(o)*C;
            Yc=2i*omega(o)*C;
            Th=[diag([cos(k*l1/2),cos(k*l2/2),cos(k*l1/2)]),1i*diag([Z01,Z02,Z01]).*diag([sin(k*l1/2),sin(k*l2/2),sin(k*l1/2)]);
            1i*diag([sin(k*l1/2),sin(k*l2/2),sin(k*l1/2)]).*diag([1/Z01,1/Z02,1/Z01]),diag([cos(k*l1/2),cos(k*l2/2),cos(k*l1/2)])];
            X=[Ya,Yb,0;
                Yb,Yc,Yb;
                0,Yb,Ya];
            Tv=[eye(3),zeros(3);
                X,eye(3)];
            T=Th*Tv*Th;
            E=diag([x1,y,x2,x1,y,x2]);
            if abs(det(T-E)) < eps
                KX=[KX,kx(i)];
                KY=[KY,ky(j)];
                F=[F,omega(o)/2/pi];
            end
        end
    end
end

delete(gcp)


%% 绘图
h1=figure;
tri=delaunay(KX,KY);
trisurf(tri,KX,KY,F);
xlabel('k1a')
ylabel('k2a')
zlabel('omega')

h2=figure;
plot3(KX,KY,F,'b.')
xlabel('k1a')
ylabel('k2a')
zlabel('omega')