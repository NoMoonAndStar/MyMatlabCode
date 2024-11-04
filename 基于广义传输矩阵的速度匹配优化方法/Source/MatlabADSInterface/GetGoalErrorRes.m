function GoalError = GetGoalErrorRes(OptimInput, ADSInterface, DatasetFile)
    %% 根据给定参数计算距离目标的误差（只优化阻抗，不优化速度）    
    ZPort1Real = OptimInput(1);
    ZPort2Real = OptimInput(2);
    ZPort4Real = OptimInput(3);
    ZPort5Real = OptimInput(4);
    % Lseries1 = OptimInput(5);
    % Lseries2 = OptimInput(6);
    Zc1 = OptimInput(5);
    Zc2 = OptimInput(6);
    % Lseries1 = OptimInput(1);
    % Lseries2 = OptimInput(2);
    % Zc1 = OptimInput(3);
    % Zc2 = OptimInput(4);
    % 修改参数并仿真
    RunSimWithParameter(ADSInterface, 'ZPort1Real', ZPort1Real, 'ZPort2Real', ZPort2Real, 'ZPort4Real', ZPort4Real, ...
        'ZPort5Real', ZPort5Real, 'Zc1', Zc1, 'Zc2', Zc2);
    % RunSimWithParameter(ADSInterface, 'Lseries1', Lseries1, 'Lseries2', Lseries2, 'Zc1', Zc1, 'Zc2', Zc2);
    
    % ZPort1Real = OptimInput(1);
    % ZPort2Real = OptimInput(2);
    % ZPort4Real = OptimInput(3);
    % ZPort5Real = OptimInput(4);
    % Ev1 = OptimInput(5);
    % Ev2 = OptimInput(6);
    % Ev4 = OptimInput(7);
    % Ev5 = OptimInput(8);
    % RunSimWithParameter(ADSInterface, 'ZPort1Real', ZPort1Real, 'ZPort2Real', ZPort2Real, 'ZPort4Real', ZPort4Real, ...
    %     'ZPort5Real', ZPort5Real, 'Ev1', Ev1, 'Ev2', Ev2, 'Ev4', Ev4, 'Ev5', Ev5)
    % 读取结果
    ADSInterface.ReadDataset(DatasetFile);

    %% 在这里定义需要优化的目标函数

    Goal = zeros(1,4);
    % weight = [3, 1, 1, 1, 1];
    weight = [2, 2, 2, 2];
    % 读取S参数
    [S41, Freq] = ADSInterface.GetVariableAsFunction('S[4,1]', 'freq');
    [S52, ~] = ADSInterface.GetVariableAsFunction('S[5,2]', 'freq');
    [S44, ~] = ADSInterface.GetVariableAsFunction('S[4,4]', 'freq');
    [S55, ~] = ADSInterface.GetVariableAsFunction('S[5,5]', 'freq');

    ftr = 10;
    fm = 200;
    % |S44| <= 0.1;
    if abs(S44(fm+1)) <= 0.1
        Goal(1) = 0;
    else
        Goal(1) = (abs(S44(fm+1))-0.1)/0.9;
    end

    % |S55| <= 0.1;
    if abs(S55(ftr+1)) <= 0.1
        Goal(2) = 0;
    else
        Goal(2) = (abs(S55(ftr+1))-0.1)/0.9;
    end

    % |S41| >= 0.99;
    if abs(S41(fm+1)) > 0.99
        Goal(3) = 0;
    else
        Goal(3) = (0.99-abs(S41(fm+1)))/0.99;
    end

    % |S52| >= 0.99;
    if abs(S52(ftr+1)) > 0.99
        Goal(4) = 0;
    else
        Goal(4) = (0.99-abs(S52(ftr+1)))/0.99;
    end
    GoalError = sum(Goal.*weight);
end