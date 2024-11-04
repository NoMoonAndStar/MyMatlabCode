%% Matlab-ADS-Interface 控制脚本
clc
clear
close all

%% 初始化

% 重要！！！必须把包含接口类的文件夹添加到路径！！！
addpath D:\Documents\GitHub\MyMatlabCode\基于广义传输矩阵的速度匹配优化方法\Source\TADSInterface

% 定义ADS接口对象
ADS = TADSInterface();
% 设置ADS安装目录
ADSInstallationDir = 'D:\Program Files\Keysight\ADS2022';
ADS.SetADSPaths(ADSInstallationDir);

% 设置参数文件的路径（以netlist.log命名）
ADS.NetlistFile = 'C:\Users\35406\Downloads\ADS-Matlab-Interface-master\ADS-Matlab-Interface-master\Demos\ADSProjects\Test_wrk\netlist.log';

% 设置仿真结果文件的路径（注意！！！直接用IDE仿真的结果和用Matlab接口仿真的结果文件所在不同！！！）
DatasetFile = 'C:\Users\35406\Downloads\ADS-Matlab-Interface-master\ADS-Matlab-Interface-master\Demos\ADSProjects\Test_wrk\test.ds';

ADS.MessageLevelOfDetails = 2;
%% 修改参数并仿真分析结果

% 修改参数并仿真
RunSimWithParameter(ADS, 'ZPort1Real', 5, 'ZPort2Real', 50);

% 读取仿真结果
ADS.ReadDataset(DatasetFile);

% 读取S21，作为频率的函数
[S41, Freq] = ADS.GetVariableAsFunction('S[4,1]', 'freq');
figure;
plot(Freq,20*log10(abs(S41)))

figure;
plot(Freq,angle(S41))
