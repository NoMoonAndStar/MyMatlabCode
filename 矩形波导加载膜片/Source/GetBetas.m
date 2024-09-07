%% 根据频率给出空间频率
function Beta = GetBetas(f)

    %初始化
    a = 60 * 1e-3; %矩形波导的长度
    b = 5.08 * 1e-3; %矩形波导的宽度
    h = 3.81 * 1e-3; %膜片高度
    p = 10 * 1e-3; %膜片间距
    omega = 2 * pi * f * 1e9;
    c = 3e8; %真空光速
    k0 = omega / c; %真空波数
    Beta10 = sqrt(k0 .^ 2 - (pi / a) ^ 2); %TE10模相位常数
    B = 4 * b * Beta10 * log(sec((pi * h) / (2 * b))) / pi; %膜片的归一化导纳

    if imag(Beta10) ~= 0
        Beta = [];
    else
        solution0 = acos(cos(Beta10 * p) - B * sin(Beta10 * p) / 2);

        if isreal(solution0)
            solutions = [];

            for j = 0:3 % 因为周期是2π，所以4π内有2个完整的周期
                solutions = [solutions solution0 + 2 * pi * j]; %#ok<AGROW>
                solutions = [solutions -solution0 + 2 * pi * (j + 1)]; %#ok<AGROW>
            end

            ValidSolutions = solutions(-4 * pi <= solutions & solutions <= 4 * pi);
            Beta = ValidSolutions;
        end

    end
