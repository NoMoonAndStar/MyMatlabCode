%% reproduce_of_PRL Broadband Nonreciprocal Amplification in Luminal Metamaterials supplementary material

% clc;
% clear;
% close all;

ordern = 1;
N = ordern * 2 + 1; %展开阶数，最终的维度阶数乘二+1
n = -ordern:1:ordern; %创建后续所需要矩阵，所用到的数组
epsilon0 = 8.85e-12;
c = 3e8;
alpha = 0.04;

g = 2 * pi;
Omega = g * c;
miu = 4 * pi * 1e-7;
epsilonr1 = 1;
d1 = 0; %介质板长度
d2 = 10;
d3 = 20;
d4 = 30;

% d1=0;%介质板长度
% d2=2*pi*c/Omega;
% d3=2*pi*c/Omega*1.5;
% d4=2*pi*c/Omega*4/pi;

xunhuan = 1;
M_combined = zeros(2 * N, 2 * N);
k1 = zeros(2 * N, 2000);

%% 求解色散曲线
for omega = linspace(-2 * Omega, 2 * Omega, 2000)
    Mee = diag(-n * g);
    Mhh = diag(-n * g);
    Meh = diag(-miu * (omega + n * Omega));
    Mhe = diag(-epsilon0 * (omega + n * Omega));

    for i = -ordern:ordern

        for j = -ordern:ordern

            if i - j == 1 %对角线下面的斜对角线
                Mhe(i + ordern + 1, j + ordern + 1) = -alpha * epsilon0 * (omega + (i) * Omega);
            end

            if i - j == -1 %对角线上面的斜对角线
                Mhe(i + ordern + 1, j + ordern + 1) = -alpha * epsilon0 * (omega + (i) * Omega);
            end

        end

    end

    M_combined = [Mee, Meh; Mhe, Mhh];
    k1(:, xunhuan) = eig(M_combined);
    xunhuan = xunhuan + 1;
end

clear i;
clear j;

k1 = real(k1);
k1 = sort(k1);
figure(1)

for i = 1:2 * N
    % plot(k(i,:)*Lambda/(2*pi),linspace(-2,2,1000),'b','linewidth',2);
    plot(k1(i, :), linspace(-2 * Omega, 2 * Omega, 2000) / Omega, 'b', 'linewidth', 1);
    % plot(k(i,:)*c/(OMEGA*sqrt(5.25)),linspace(-2*c,2*c,1000)/OMEGA,'b','linewidth',1);
    hold on
    grid on
end

clear i;
clear j;
%% 构造最后的线性方程组，这里需要构造Mvinc和Mvref
% alpha=0;
omega = 0.69 * Omega;
% omega=10e6;
kv = zeros(1, 2 * N);
Hv = zeros(1, 2 * N);
kv(1, 1:N) = -n * g - sqrt(((omega + n * Omega) / c) .^ 2); %先对负号的kv进行赋值
kv(1, N + 1:2 * N) = -n * g + sqrt(((omega + n * Omega) / c) .^ 2); %再对正号的kv进行赋值

Hv(1, 1:N) =- (n * g + kv(1, 1:N)) ./ (miu * (omega + n * Omega)); %先计算前三个特征值kv所对应的磁场的特征向量的值
Hv(1, N + 1:2 * N) =- (n * g + kv(1, N + 1:2 * N)) ./ (miu * (omega + n * Omega)); %再计算后三个特征值kv所对应的磁场的特征向量的值

% Hv(1,1:N)=-1/(miu*c);%先计算前三个特征值kv所对应的磁场的特征向量的值
% Hv(1,N+1:2*N)=1/(miu*c);%再计算后三个特征值kv所对应的磁场的特征向量的值
M11 = diag(ones(1, N));
M12 = diag(ones(1, N));
M21 = diag(Hv(1, 1:N));
M22 = diag(Hv(1, N + 1:2 * N));
Mv = [M11, M12; M21, M22];
index = find(Hv < 0);
Mvinc = Mv(:, index);
Mvref = Mv;
Mvref(:, index) = [];
% Mvinc=zeros(2*N,N);
% for i=1:N
%     Mvinc(i,i)=1;
%     Mvinc(i+N,i)=-1/(miu*c);
% end
% Mvref=zeros(2*N,N);
% for i=1:N
%     Mvref(i,i)=1;
%     Mvref(i+N,i)=+1/(miu*c);
% end
% new_order=[2,3,1];
% Mvref=Mvref(:,new_order);
% eigenzheng=zeros(N,2000);
% eigenfu=zeros(N,2000);
% P=zeros(2*N,2*N);%介质传播的相移矩阵
evinc = zeros(N, 1);
evinc(ordern + 1, 1) = 1; %因为Evinc已知，入射波未平面波展开，所以在n=0时有值，这里值可以设置
% y=1;
% evtra=zeros(N,2000);
% for omega=linspace(-2*Omega,2*Omega,2000)
%     Mee=diag(-n*g);
%     Mhh=diag(-n*g);
%     Meh=diag(-miu*(omega+n*Omega));
%     Mhe=diag(-epsilon0*(omega+n*Omega));
%     for i=1:N
%         for j=1:N
%             if i-j==1%对角线下面的斜对角线
%                 Mhe(i,j)=-alpha*epsilon0*(omega+(i-2*ordern)*Omega);
%             end
%             if i-j==-1%对角线上面的斜对角线
%                 Mhe(i,j)=-alpha*epsilon0*(omega+(i-2*ordern)*Omega);
%             end
%         end
%     end
%     M_combined=[Mee,Meh;Mhe,Mhh];
%     clear i;
%     clear j;
%     [V,D] = eig(M_combined);                        %求Vaccum中的本征场，根据推导特征向量应该只有位置的不同，值是相同的
%     k1=eig(M_combined);
%     Mm=V;
%     P=diag(exp(1i*d*k1));                          %介质传播的相移矩阵
%     A=Mm*(eye(2*N,2*N)/(Mm*P))*Mvinc;
%     % A=Mm*inv(Mm*P)*Mvinc;
%     B=-Mvref;
%     coffient=Mvinc*evinc;
%     AB=[A,B];
%     ev1=AB\coffient;
%     evtra(:,y)=ev1(1:N,:);
%     eigenzheng(:,y)=-n*g+sqrt((c^-2)*(omega+n*Omega).^2);
%     eigenfu(:,y)=-n*g-sqrt((c^-2)*(omega+n*Omega).^2);
%
%     y=y+1;
% end
% clear i;
% clear j;

%% 中间段 用于计算单个频率下的透射场，和老师所选的频率进行验证。
% omega=0.69*Omega;
% omega=10e6;
Mee = diag(-n * g);
Mhh = diag(-n * g);
Meh = diag(-miu * (omega + n * Omega));
Mhe = diag(-epsilon0 * (omega + n * Omega));

for i = -ordern:ordern

    for j = -ordern:ordern

        if i - j == 1 %对角线下面的斜对角线
            Mhe(i + ordern + 1, j + ordern + 1) = -alpha * epsilon0 * (omega + (i) * Omega);
        end

        if i - j == -1 %对角线上面的斜对角线
            Mhe(i + ordern + 1, j + ordern + 1) = -alpha * epsilon0 * (omega + (i) * Omega);
        end

    end

end

M_combined = [Mee, Meh; Mhe, Mhh];
clear i;clear j;
[V, D] = eig (M_combined); %求介质板中的本征场，根据推导特征向量应该只有位置的不同，值是相同的
k1 = diag(D);
re_k1 = real(k1);
[Y, I] = sort(re_k1);
k1 = k1(I);
V = V(:, I);
Mm = V;
% Mm=Mm1;
clear i;clear j;
P1 = diag(exp(1i * d1 * k1)); %介质传播的相移矩阵
P2 = diag(exp(1i * k1 * d2));
P3 = diag(exp(1i * d3 * k1));
P4 = diag(exp(1i * d4 * k1));
A1 = Mm * ((Mm * P1) \ Mvinc);
B = -Mvref;
coffient = Mvinc * evinc;
AB1 = [A1, B];
ev1 = AB1 \ coffient;
evtra1 = ev1(1:N, :);
Etn1 = Mvinc * evtra1;
k1 = omega / c;
omegat = 0:pi / 360:4 * pi;

Etz1 = zeros(1, length(omegat));

for l = 1:N
    Etz1 = Etn1(l, 1) * exp(1i * ((k1 + (l - ordern - 1) * g) * d1 - (omega / Omega + (l - ordern - 1)) * omegat)) + Etz1;
end

% Etz1=Etn1(1,1)*exp(1i*((k-1*g)*d1-(omega/Omega-1)*omegat))+Etn1(2,1)*exp(1i*((k-0*g)*d1-(omega/Omega-0)*omegat))+Etn1(3,1)*exp(1i*((k+g)*d1-(omega/Omega-1)*omegat));%组装Etz
power1 = abs(Etz1) .^ 2;

A2 = Mm * ((Mm * P2) \ Mvinc);
coffient = Mvinc * evinc;
AB2 = [A2, B];
ev2 = AB2 \ coffient;
evtra2 = ev2(1:N, :);
Etn2 = Mvinc * evtra2;
Etz2 = zeros(1, length(omegat));
Etz2_matrix = zeros(N, length(omegat));

for l = 1:N
    Etz2_matrix(l, :) = real(Etn2(l, 1) * exp(1i * (((k1 + (l - ordern - 1) * g) * d2 - (omega / Omega + (l - ordern - 1)) * omegat))));
    Etz2 = Etn2(l, 1) * exp(1i * ((k1 + (l - ordern - 1) * g) * d2 - (omega / Omega + (l - ordern - 1)) * omegat)) + Etz2;
end

% Etz2=sum(Etn2(1:N))
% figure(6)
% plot(omegat/Omega,Etz2_matrix(1,:))
% hold on;
% plot(omegat/Omega,Etz2_matrix(2,:))
% plot(omegat/Omega,Etz2_matrix(3,:))
% plot(omegat/Omega,Etz2_matrix(4,:))
% plot(omegat/Omega,Etz2_matrix(5,:))
% plot(omegat/Omega,Etz2_matrix(6,:))
% plot(omegat/Omega,Etz2_matrix(7,:))
% plot(omegat/Omega,Etz2_matrix(8,:))
% plot(omegat/Omega,Etz2_matrix(9,:))
% Etz2=Etn2(1,1)*exp(1i*((k-1*g)*d2-(omega/Omega-1)*omegat))+Etn2(2,1)*exp(1i*((k-0*g)*d2-(omega/Omega-0)*omegat))+Etn2(3,1)*exp(1i*((k+g)*d2-(omega/Omega+1)*omegat));%组装Etz
power2 = abs(Etz2) .^ 2;

A3 = Mm * (eye(2 * N, 2 * N) / (Mm * P3)) * Mvinc;
coffient = Mvinc * evinc;
AB3 = [A3, B];
ev3 = AB3 \ coffient;
evtra3 = ev3(1:N, :);
Etn3 = Mvinc * evtra3;

Etz3 = zeros(1, length(omegat));

for l = 1:N
    Etz3 = Etn3(l, 1) * exp(1i * ((k1 + (l - ordern - 1) * g) * d3 - (omega / Omega + (l - ordern - 1)) * omegat)) + Etz3;
end

% for l=1:N
%     Etz3=Etn3(l,1)*exp(1i*((k)*d3-(omega/Omega+(l-ordern-1))*omegat))+Etz3;
% end

% Etz3=Etn3(1,1)*exp(1i*((k-1*g)*d3-(omega/Omega-1)*omegat))+Etn3(2,1)*exp(1i*((k-0*g)*d3-(omega/Omega-0)*omegat))+Etn3(3,1)*exp(1i*((k+g)*d3-(omega/Omega+1)*omegat));%组装Etz
power3 = abs(Etz3) .^ 2;

% A4=Mm*(eye(2*N,2*N)/(Mm*P4))*Mvinc;
A4 = Mm * ((Mm * P4) \ Mvinc);
coffient = Mvinc * evinc;
AB4 = [A4, B];
ev4 = AB4 \ coffient;
evtra4 = ev4(1:N, :);
Etn4 = Mvinc * evtra4;
Etz4 = zeros(1, length(omegat));

for l = 1:N
    Etz4 = Etz4 + Etn4(l, 1) * exp(1i * (((k1 + (l - ordern - 1) * g) * d4 - (omega / Omega + (l - ordern - 1)) * omegat)));
end

% Etz4=Etn4(1,1)*exp(1i*((k-1*g)*d4-(omega/Omega-1)*omegat))+Etn4(2,1)*exp(1i*((k-0*g)*d4-(omega/Omega-0)*omegat))+Etn4(3,1)*exp(1i*((k+g)*d4-(omega/Omega+1)*omegat));%组装Etz
power4 = abs(Etz4) .^ 2;
figure(2)
% subplot(2,1,1)
plot(omegat, power1 / max(power3), 'r'); %max(power4)
% plot(omegat,power1,'r');
hold on;
plot(omegat, power2 / max(power3), 'b');
plot(omegat, power3 / max(power3), 'g');
plot(omegat, power4 / max(power3), 'k');
legend('|E_{out}|^2(d=0) ', '|E_{out}|^2(d=10) ', '|E_{out}|^2(d=20) ', '|E_{out}|^2(d=30) ');
title('Transmitted wave')
ylabel('|E_{out}|^2')
xlabel('\Omegat')
% subplot(2,1,2)
figure(3)
epsilonr1 = 1 + 2 * alpha * cos(g * d1 - omegat);
hold on;
plot(omegat, epsilonr1);
epsilonr2 = 1 + 2 * alpha * cos(g * d2 - omegat);
plot(omegat, epsilonr2);
epsilonr3 = 1 + 2 * alpha * cos(g * d3 - omegat);
plot(omegat, epsilonr3);
epsilonr4 = 1 + 2 * alpha * cos(g * d4 - omegat);
plot(omegat, epsilonr4);
legend('d=0 ', 'd=10 ', 'd=20 ', 'd=30 ');
ylabel('epslionr')
xlabel('\Omegat')

%% 分界线
% Etn=Mvinc*evtra;% 求出透射电场的系数，然后对其选取频率，在对展开的阶数进行叠加
% %% 选取频率，通过叠加解出透射场
% omega=linspace(-2*Omega,2*Omega,2000);
% index=find(omega==0.69*Omega);
% omega(1351)% 选第700个点
% k=omega(700)/c;
% omegat=0:0.001:6.5*pi;
% t=omegat/Omega;
% Etz=Etn(1,1351)*exp(1i*((k-1*g)*d-(omega(1351)/Omega-1)*omegat))+Etn(2,1351)*exp(1i*((k-0*g)*d-(omega(1351)/Omega-0)*omegat))+Etn(3,1351)*exp(1i*((k+g)*d-(omega(1351)/Omega-1)*omegat));%组装Etz
% fuzhi=abs(Etz);
% power=fuzhi.^2;
% figure(2)
% subplot(2,1,1)
% plot(omegat,power/max(power));
% % xline(0.5*pi, 'r--')
% % % xticks(0.5*pi); % 设置横坐标刻度为 0.5pi
% % % xticklabels({'0.5\pi'}); % 设置横坐标标签为 '0.5pi'】
% % xline(1.5*pi, 'r--')
% subplot(2,1,2)
% epsilont=1+2*alpha*cos(g*d-omegat);
% plot(omegat,epsilont);
% % xline(1.5*pi, 'r--')
% % xline(0.5*pi, 'r--')
%
%
%
%

% t=linspace(0,10*pi/Omega,1000);
% x=d;
% T=Omega*t;
% X=x-T/g;
% U=exp(-2*alpha*g*d*sin(g*X)-(alpha*g*d*cos(g*X)).^2);
% figure(3);
% plot(T,U);

% t=linspace(0,200e-9,1000);
% Etz4=zeros(1,length(t));
% for l=1:N
%     Etz4=Etz4+Etn4(l,1)*exp(1i*(((k+(l-ordern-1)*g)*d4-(omega+(l-ordern-1)*Omega)*t)));
% end
%
% figure(4)
% plot(t*1e9,abs(Etz4).^2,'linewidth',2);
% xlabel('ns');
% ylabel('|E^2|');

figure(8)
plot(omegat / Omega, real(Etz1), 'k', 'linewidth', 2)
hold on;
% t=omegat/3e8;
% plot(t,cos(0.69*OMEGA*t))
xlabel('t(s)')
ylabel('magnitude')

clear i;
Etz2_time = zeros(1, length(omegat));

for i = 1:N
    Etz2_time = real(Etz2_matrix(i, :)) + Etz2_time;
end

figure(9)
plot(omegat / Omega, real(Etz2_time), 'b', 'linewidth', 2)
hold on;
plot(omegat / Omega, real(Etz2), 'k', 'linewidth', 1)
xlabel('t(s)')
ylabel('magnitude')

figure(10)

plot(omegat / Omega, real(Etz3), 'r', 'linewidth', 2)
hold on;
xlabel('t(s)')
ylabel('magnitude')

figure(11)
Etz4_real = real(Etz4);
plot(omegat / Omega, real(Etz4) / max(Etz4_real), 'g', 'linewidth', 2)
hold on;
xlabel('t(s)')
ylabel('magnitude')
title('d=30')
