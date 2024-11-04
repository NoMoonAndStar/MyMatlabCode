function [ str_varName ] = getVarName( var )
%% 将变量的变量名转换为字符串

    str_varName = sprintf('%s', inputname(1));

end