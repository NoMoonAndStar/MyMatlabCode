function GoalError = TimeGetGoalError(OptimInput, ADSInterface, DatasetFile)
    %% 根据给定参数计算距离目标的误差    
    ZPort1 = OptimInput(1);
    ZPort2 = OptimInput(2);
    ZPort4 = OptimInput(3);
    ZPort5 = OptimInput(4);
    Zc1 = OptimInput(5);
    Zc2 = OptimInput(6);
    % 修改参数并仿真
    RunSimWithParameter(ADSInterface, 'ZPort1', ZPort1, 'ZPort2', ZPort2, 'ZPort4', ZPort4, ...
        'ZPort5', ZPort5, 'Zc1', Zc1, 'Zc2', Zc2);
    ADSInterface.ReadDataset(DatasetFile);

    %% 在这里定义需要优化的目标函数

    [m1, time] = ADSTime.GetVariableAsFunction('m1', 'time');
    % 主路波形
    trwave = zeros(length(time), 33);
    for i = 0:32
        [trwave(:,i+1), ~] = ADSTime.GetVariableAsFunction(strcat('tr',num2str(i)), 'time');
    end

    % 计算理论放大曲线
    Cj0 = 8.47*1e-12;
    M = 70;
    Vj0 = 80;
    V_static = 3;
    a = 5*1e-3;
    epsr = 1.3;
    m = (max(m1) + V_static) / V_static;
    fm = 200*1e6;
    
    [max_val1]

    t = 
    [~, TheoryPA] = TheoryAmplified(Cj0, M, Vj0, V_static, a, epsr, Zc2, m, fm);
end