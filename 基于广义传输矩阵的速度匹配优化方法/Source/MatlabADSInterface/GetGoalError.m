function GoalError = GetGoalError(OptimInput, ADSInterface, DatasetFile)
    %% 根据给定参数计算距离目标的误差    
    ZPort1 = OptimInput(1);
    ZPort2 = OptimInput(2);
    ZPort4 = OptimInput(3);
    ZPort5 = OptimInput(4);
    Zc1 = OptimInput(5);
    Zc2 = OptimInput(6);

    RunSimWithParameter(ADSInterface, 'ZPort1', ZPort1, 'ZPort2', ZPort2, 'ZPort4', ZPort4, ...
        'ZPort5', ZPort5, 'Zc1', Zc1, 'Zc2', Zc2);

    ADSInterface.ReadDataset(DatasetFile);

    %% 在这里定义需要优化的目标函数

    Goal = zeros(1,17);
    % weight = [3, 1, 1, 1, 1];
    weight = [2, 2, 8, 2, 1, 1, 1, 1, 1, 1 60, 1, 1, 1, 1, 1, 1];
    % weight = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 1, 1, 0.5, 0.5, 0.25, 0.25];
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

    %|S52| >= 0.99;
    if abs(S52(fm-ftr+1)) > 0.99
        Goal(5) = 0;
    else
        Goal(5) = (0.99-abs(S52(fm-ftr+1)))/0.99;
    end

    % |S52| >= 0.99;
    if abs(S52(fm+ftr+1)) > 0.99
        Goal(6) = 0;
    else
        Goal(6) = (0.99-abs(S52(fm+ftr+1)))/0.99;
    end

    % |S52| >= 0.99;
    if abs(S52(2*fm-ftr+1)) > 0.99
        Goal(7) = 0;
    else
        Goal(7) = (0.99-abs(S52(2*fm-ftr+1)))/0.99;
    end

    % |S52| >= 0.99;
    if abs(S52(2*fm+ftr+1)) > 0.99
        Goal(8) = 0;
    else
        Goal(8) = (0.99-abs(S52(2*fm+ftr+1)))/0.99;
    end
    
    % |S52| >= 0.99;
    if abs(S52(3*fm-ftr+1)) > 0.99
        Goal(9) = 0;
    else
        Goal(9) = (0.99-abs(S52(2*fm+ftr+1)))/0.99;
    end

    % |S52| >= 0.99;
    if abs(S52(3*fm+ftr+1)) > 0.99
        Goal(10) = 0;
    else
        Goal(10) = (0.99-abs(S52(2*fm+ftr+1)))/0.99;
    end

    % 速度匹配;
    delta41 = phase(S41)./Freq/2/pi;
    delta52 = phase(S52)./Freq/2/pi;
    Goal(11) = delta41(fm+1) -delta52(ftr+1);
    Goal(11) = abs(Goal(11) / (Goal(11) + delta41(fm+1)));
    Goal(12) = delta41(fm+1) - delta52(fm+ftr+1);    % 一阶
    Goal(12) = abs(Goal(12) / (Goal(12) + delta41(fm+1)));
    Goal(13) = delta41(fm+1) - delta52(fm-ftr+1);    % 一阶
    Goal(13) = abs(Goal(13) / (Goal(13) + delta41(fm+1)));
    Goal(14) = delta41(fm+1) - delta52(2*fm-ftr+1);  % 二阶
    Goal(14) = abs(Goal(14) / (Goal(14) + delta41(fm+1)));
    Goal(15) = delta41(fm+1) - delta52(2*fm+ftr+1);  % 二阶
    Goal(15) = abs(Goal(15) / (Goal(15) + delta41(fm+1)));
    Goal(16) = delta41(fm+1) - delta52(3*fm-ftr+1);  % 三阶
    Goal(16) = abs(Goal(16) / (Goal(16) + delta41(fm+1)));
    Goal(17) = delta41(fm+1) - delta52(3*fm+ftr+1);  % 三阶
    Goal(17) = abs(Goal(17) / (Goal(17) + delta41(fm+1)));
    GoalError = sum(Goal.*weight);
end