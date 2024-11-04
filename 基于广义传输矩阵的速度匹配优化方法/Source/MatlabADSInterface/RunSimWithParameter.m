function [] = RunSimWithParameter(ADSInterface, varargin)
%修改ADS的参数，并运行;
%
%   ADSInterface：ADS接口，必须是TADSInterface的实例，
%   需要定义好安装目录、数据文件位置、参数文件位置；
%
%   varargin：参数名，参数值对，参数名必须为字符串类型，
%   参数值必须为数值类型；

    % 读取参数
    ParameterName = varargin(1:2:end);  % 要操作的参数名称
    ParameterValue = varargin(2:2:end); % 要操作的参数值

    % 参数验证
    mustBeA(ADSInterface, 'TADSInterface') % ADS与matlab接口对象

    if (isempty(ParameterValue))
    elseif ~(length(ParameterName)==length(ParameterValue))
        error('输入参数名与参数值不匹配')
    elseif ~(ischar(cell2mat(ParameterName)))
        error('参数名必须为字符串类型')
    elseif ~(isnumeric(cell2mat(ParameterValue)))
        error('参数值必须为数值类型')
    end

    % 修改ADS参数
    for ip = 1:length(ParameterName)
        ADSInterface.ChangeParameter(cell2mat(ParameterName(ip)), cell2mat(ParameterValue(ip)),'double');
    end

    % 运行仿真
    ADSInterface.RunSimulation();
end