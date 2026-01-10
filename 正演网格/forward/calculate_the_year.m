clc; clear;

%% 1. 参数定义
d13C_volcanic = -26.0;   % 火山 CO2 同位素
d13C_methane  = -35.0;  % 甲烷同位素
d13C_target   = -27.8;  % 原本设定的混合目标值

%% 2. 数据加载
pars.time = [ ...
    -3.00e+06 0.00e+00 0.1 5.00e+01 1.50e+02 2.50e+02 3.50e+02 4.50e+02 5.50e+02 ...
    6.50e+02 7.50e+02 8.50e+02 9.50e+02 1.05e+03 1.15e+03 1.25e+03 1.35e+03 ...
    1.45e+03 1.55e+03 1.65e+03 1.75e+03 1.85e+03 1.95e+03 2.05e+03 2.15e+03 ...
    2.25e+03 2.35e+03 2.45e+03 2.55e+03 2.65e+03 2.75e+03 2.85e+03 2.95e+03 ...
    3.05e+03 3.15e+03 3.25e+03 3.35e+03 3.45e+03 3.55e+03 3.65e+03 3.75e+03 ...
    3.85e+03 3.95e+03 4.05e+03 4.15e+03 4.25e+03 4.35e+03 4.45e+03 ...
    4.55e+03 4.65e+03 4.75e+03 4.85e+03 4.95e+03 5.05e+03 5.15e+03 5.25e+03 ...
    5.35e+03 5.45e+03 5.55e+03 5.65e+03 5.75e+03 5.85e+03 5.95e+03 6000];

pars.input_C = [ ...
    0.00e+00 0.00e+00 2.52e13 2.72e+13 3.12e+13 3.42e+13 3.66e+13 3.81e+13 4.00e+13 ...
    4.16e+13 4.31e+13 4.45e+13 4.55e+13 4.69e+13 4.77e+13 4.88e+13 4.99e+13 ...
    5.06e+13 5.16e+13 5.25e+13 5.36e+13 5.46e+13 5.57e+13 5.66e+13 5.78e+13 ...
    5.87e+13 5.97e+13 6.06e+13 6.19e+13 6.31e+13 6.44e+13 6.55e+13 6.67e+13 ...
    6.78e+13 6.93e+13 7.03e+13 7.14e+13 7.29e+13 7.39e+13 7.48e+13 7.58e+13 ...
    7.82e+13 7.97e+13 8.11e+13 8.24e+13 8.37e+13 8.52e+13 8.68e+13 8.83e+13 ...
    8.98e+13 9.13e+13 9.31e+13 9.43e+13 9.65e+13 9.83e+13 1.00e+14 1.02e+14 ...
    1.04e+14 1.06e+14 1.08e+14 1.10e+14 1.12e+14 1.13e+14 1.135e14 ] ;

total_len = length(pars.time); % 获取总长度 (例如 62)

%% 3. 分割设置
split_index = 58; % 你选择的分割点 (第54个点开始变为纯火山)

% 验证索引有效性
if split_index > total_len || split_index < 3
    error('分割索引超出范围');
end
fprintf('分割点: Index %d, 时间 %.1f yr\n', split_index, pars.time(split_index));

%% 4. 积分计算总碳量 (使用 trapz 自动处理非等距步长)
% 提取 Spin-up 之后的数据 (从第3个点开始)
valid_indices = 3:total_len;
time_event = pars.time(valid_indices);
rate_event = pars.input_C(valid_indices);

% trapz(time, rate) = 积分 (Rate * dt)，这才是真正的物质总量
Total_Carbon_Moles = trapz(time_event, rate_event);

%% 5. 反推原本的物质分配 (Target: -27.8)
% 杠杆原理: -27.8 = f_CH4 * (-35) + (1-f_CH4) * (-6)
frac_CH4_global = (d13C_target - d13C_volcanic) / (d13C_methane - d13C_volcanic);

Moles_CH4_Total = Total_Carbon_Moles * frac_CH4_global;
Moles_CO2_Total = Total_Carbon_Moles * (1 - frac_CH4_global);

fprintf('=== 原本预算 ===\n');
fprintf('总排放量: %.4e mol\n', Total_Carbon_Moles);
fprintf('甲烷总量: %.4e mol (占比 %.2f%%)\n', Moles_CH4_Total, frac_CH4_global*100);
fprintf('CO2总量:  %.4e mol\n', Moles_CO2_Total);

%% 6. 分段积分计算
% 阶段1: 混合段 (indices 3 到 split_index-1)
% 阶段2: 纯火山段 (indices split_index 到 end)

% 注意积分区间的衔接：
% 严格来说 trapz(3:end) = trapz(3:split) + trapz(split:end) 
% 这里的分割点属于两个区间，但在模型赋值时，该点的值属于后段
idx_p1 = 3 : split_index;      % 为了积分准确，包含边界点
idx_p2 = split_index : total_len; 

% 计算各段的物理总量 (Mass)
Mass_P1_Raw = trapz(pars.time(idx_p1), pars.input_C(idx_p1));
Mass_P2_Raw = trapz(pars.time(idx_p2), pars.input_C(idx_p2));

% 稍微修正：由于离散积分边界重复计算了一次边界点的一半，
% 我们通过比例归一化以确保 Mass_P1 + Mass_P2 正好等于 Total
scale_factor = Total_Carbon_Moles / (Mass_P1_Raw + Mass_P2_Raw);
Mass_P1 = Mass_P1_Raw * scale_factor;
Mass_P2 = Mass_P2_Raw * scale_factor;

fprintf('\n=== 分段分配 ===\n');
fprintf('后段 (纯火山) 预计排放量: %.4e mol\n', Mass_P2);

%% 7. 物质再分配
% 后段全是 CO2
Moles_CH4_P2 = 0;
Moles_CO2_P2 = Mass_P2;

% 检查是否还有足够的 CO2 给后段
if Moles_CO2_P2 > Moles_CO2_Total
    error('错误：后段所需的碳量超过了原本模型中火山CO2的总预算！');
end

% 前段承担所有甲烷 + 剩余 CO2
Moles_CH4_P1 = Moles_CH4_Total; 
Moles_CO2_P1 = Moles_CO2_Total - Moles_CO2_P2;
Mass_P1_Calc = Moles_CH4_P1 + Moles_CO2_P1;

%% 8. 计算前段新同位素值
d13C_New_Phase1 = (Moles_CH4_P1 * d13C_methane + Moles_CO2_P1 * d13C_volcanic) / Mass_P1_Calc;

fprintf('前段 (混合段) 总量: %.4e mol\n', Mass_P1_Calc);
fprintf('前段所需新同位素值: %.4f\n', d13C_New_Phase1);

%% 9. 自动生成代码 (动态长度)
fprintf('\n============================================\n');
fprintf('请复制以下代码到你的脚本中 (已自动适配长度):\n');
fprintf('============================================\n\n');

% 计算向量各部分长度
len_spinup = 2;                     % -3e6 和 0
len_p1     = split_index - 1 - 2;   % 混合段长度 (从 index 3 到 split-1)
len_p2     = total_len - split_index + 1; % 火山段长度 (从 split 到 end)

fprintf('%% === Modified Isotope Input (Auto-Calculated) ===\n');
fprintf('mu_mixed  = %.4f; %% 计算所得前段值\n', d13C_New_Phase1);
fprintf('mu_volc   = %.1f;   %% 后段纯火山值\n', d13C_volcanic);
fprintf('spin_val  = -25;    %% Spin-up 值\n\n');

fprintf('%% 自动构建向量 (总长度 %d)\n', total_len);
fprintf('pars.d13C_input = [spin_val, spin_val, ...\n');
fprintf('                   ones(1, %d) * mu_mixed, ...\n', len_p1);
fprintf('                   ones(1, %d) * mu_volc];\n\n', len_p2);

fprintf('%% 检查长度 (可选)\n');
fprintf('if length(pars.d13C_input) ~= length(pars.time)\n');
fprintf('    error(''长度不匹配! Input:%%d, Time:%%d'', length(pars.d13C_input), length(pars.time));\n');
fprintf('end\n');