%% 计算加权 NLL 指标 (平衡 准度、精度 与 k权重)
% 修正版：分别寻找每个文件夹内 p最优、Delta最优、Mu最优 的文件

clear; clc;

%% ========== 1. 参数设置 ==========
k_weight = 0.2; 

%% ========== 2. 读入 reference 数据 ==========
ref_d13C = load('data/d13Cp_data.txt');
ref_CO2  = load('data/CO2.txt');
ref_SST  = load('data/SST.txt');

ref.d13C = make_ref_struct(ref_d13C);
ref.CO2  = make_ref_struct(ref_CO2);
ref.SST  = make_ref_struct(ref_SST);

%% ========== 3. 遍历工作目录 ==========
rootDir = 'D:\desktop\waiting';
dirInfo = dir(rootDir);

% 初始化三个独立的结果容器
summary_p = struct('folder', {}, 'param1', {}, 'param2', {}, 'csv_name', {}, 'c', {}, 'd', {}, 'score_val', {}, 'other_metric1', {}, 'other_metric2', {});
summary_delta = summary_p;
summary_mu = summary_p;

for iDir = 1:numel(dirInfo)
    if ~dirInfo(iDir).isdir, continue; end
    folderName = dirInfo(iDir).name;
    if startsWith(folderName, '.') || ~contains(folderName, '_'), continue; end

    folderPath = fullfile(rootDir, folderName);
    fprintf('\n=== 处理文件夹: %s ===\n', folderName);

    % 解析文件夹参数
    folderTokens = regexp(folderName, '(-?[\d\.]+)_(-?[\d\.]+)', 'tokens', 'once');
    if ~isempty(folderTokens)
        f_p1 = str2double(folderTokens{1});
        f_p2 = str2double(folderTokens{2});
    else
        f_p1 = NaN; f_p2 = NaN;
    end

    csvList = dir(fullfile(folderPath, 'CHEES_*.csv'));
    if isempty(csvList), continue; end

    % --- 初始化三个维度的最佳记录 ---
    best_val_p = Inf;     rec_p = [];
    best_val_delta = Inf; rec_delta = [];
    best_val_mu = Inf;    rec_mu = [];

    for iFile = 1:numel(csvList)
        csvName = csvList(iFile).name;
        csvPath = fullfile(folderPath, csvName);
        
        if csvList(iFile).bytes < 10*1024 || csvList(iFile).bytes > 500*1024
            continue; 
        end

        try
            T = readtable(csvPath);
            if height(T) < 2, continue; end

            i0 = 1; iT = height(T);

            recon.d13C = make_recon_struct(T, i0, iT, 'd13C_DIC_mean', 'd13C_DIC_p5', 'd13C_DIC_p95');
            recon.CO2  = make_recon_struct(T, i0, iT, 'CO2_mean_ppm', 'CO2_p5_ppm', 'CO2_p95_ppm');
            recon.SST  = make_recon_struct(T, i0, iT, 'SST_mean_degC', 'SST_p5_degC', 'SST_p95_degC');

            % 计算 NLL 分量
            [sDelta_1, sMu_1] = compute_nll_components(ref.d13C, recon.d13C);
            [sDelta_2, sMu_2] = compute_nll_components(ref.CO2,  recon.CO2);
            [sDelta_3, sMu_3] = compute_nll_components(ref.SST,  recon.SST);

            % 平均分
            curr_delta = mean([sDelta_1, sDelta_2, sDelta_3], 'omitnan');
            curr_mu    = mean([sMu_1,    sMu_2,    sMu_3],    'omitnan');

            if isnan(curr_delta) || isnan(curr_mu), continue; end

            % 综合分
            curr_p = curr_delta + (1 / k_weight^2) * curr_mu;
            %curr_p = curr_delta * curr_mu;

            % 解析 c, d
            numStrs = regexp(csvName, '-?[\d\.]+', 'match');
            if numel(numStrs) >= 2
                curr_c = str2double(numStrs{end-1});
                curr_d = str2double(numStrs{end});
            else
                curr_c = NaN; curr_d = NaN;
            end

            % --- [核心修改] 独立更新三个最佳记录 ---
            
            % 1. 更新综合 P 最优
            if curr_p < best_val_p
                best_val_p = curr_p;
                rec_p.csv = csvName; rec_p.c = curr_c; rec_p.d = curr_d;
                rec_p.score = curr_p;
                rec_p.delta = curr_delta; rec_p.mu = curr_mu; % 附带记录其他指标
            end

            % 2. 更新 Delta 最优
            if curr_delta < best_val_delta
                best_val_delta = curr_delta;
                rec_delta.csv = csvName; rec_delta.c = curr_c; rec_delta.d = curr_d;
                rec_delta.score = curr_delta;
                rec_delta.p = curr_p; rec_delta.mu = curr_mu;
            end

            % 3. 更新 Mu 最优
            if curr_mu < best_val_mu
                best_val_mu = curr_mu;
                rec_mu.csv = csvName; rec_mu.c = curr_c; rec_mu.d = curr_d;
                rec_mu.score = curr_mu;
                rec_mu.p = curr_p; rec_mu.delta = curr_delta;
            end

        catch
            continue;
        end
    end

    % --- 保存结果到对应的 Summary ---
    if ~isempty(rec_p)
        idx = length(summary_p) + 1;
        summary_p(idx).folder = folderName;
        summary_p(idx).param1 = f_p1; summary_p(idx).param2 = f_p2;
        summary_p(idx).csv_name = rec_p.csv;
        summary_p(idx).c = rec_p.c; summary_p(idx).d = rec_p.d;
        summary_p(idx).score_val = rec_p.score;       % 这里的 score 是综合 p
        summary_p(idx).other_metric1 = rec_p.delta;   % 顺便记下它的 delta
        summary_p(idx).other_metric2 = rec_p.mu;      % 顺便记下它的 mu
    end

    if ~isempty(rec_delta)
        idx = length(summary_delta) + 1;
        summary_delta(idx).folder = folderName;
        summary_delta(idx).param1 = f_p1; summary_delta(idx).param2 = f_p2;
        summary_delta(idx).csv_name = rec_delta.csv;
        summary_delta(idx).c = rec_delta.c; summary_delta(idx).d = rec_delta.d;
        summary_delta(idx).score_val = rec_delta.score; % 这里的 score 是 delta
        summary_delta(idx).other_metric1 = rec_delta.p; 
        summary_delta(idx).other_metric2 = rec_delta.mu;
    end

    if ~isempty(rec_mu)
        idx = length(summary_mu) + 1;
        summary_mu(idx).folder = folderName;
        summary_mu(idx).param1 = f_p1; summary_mu(idx).param2 = f_p2;
        summary_mu(idx).csv_name = rec_mu.csv;
        summary_mu(idx).c = rec_mu.c; summary_mu(idx).d = rec_mu.d;
        summary_mu(idx).score_val = rec_mu.score;     % 这里的 score 是 mu
        summary_mu(idx).other_metric1 = rec_mu.p;
        summary_mu(idx).other_metric2 = rec_mu.delta;
    end
end

%% ========== 4. 排序并输出 CSV ==========

% 1. P Summary (按 P 值升序)
if ~isempty(summary_p)
    T_p = struct2table(summary_p);
    T_p = sortrows(T_p, 'score_val', 'ascend');
    % 重命名列以便阅读
    T_p.Properties.VariableNames{'score_val'} = 'Best_P';
    T_p.Properties.VariableNames{'other_metric1'} = 'Delta_at_BestP';
    T_p.Properties.VariableNames{'other_metric2'} = 'Mu_at_BestP';
    writetable(T_p, 'summary_weighted_p.csv');
    fprintf('已生成 summary_weighted_p.csv\n');
end

% 2. Delta Summary (按 Delta 值升序)
if ~isempty(summary_delta)
    T_d = struct2table(summary_delta);
    T_d = sortrows(T_d, 'score_val', 'ascend');
    T_d.Properties.VariableNames{'score_val'} = 'Best_Delta';
    T_d.Properties.VariableNames{'other_metric1'} = 'P_at_BestDelta';
    T_d.Properties.VariableNames{'other_metric2'} = 'Mu_at_BestDelta';
    writetable(T_d, 'summary_delta_only.csv');
    fprintf('已生成 summary_delta_only.csv\n');
end

% 3. Mu Summary (按 Mu 值升序)
if ~isempty(summary_mu)
    T_m = struct2table(summary_mu);
    T_m = sortrows(T_m, 'score_val', 'ascend');
    T_m.Properties.VariableNames{'score_val'} = 'Best_Mu';
    T_m.Properties.VariableNames{'other_metric1'} = 'P_at_BestMu';
    T_m.Properties.VariableNames{'other_metric2'} = 'Delta_at_BestMu';
    writetable(T_m, 'summary_mu_only.csv');
    fprintf('已生成 summary_mu_only.csv\n');
end

% ... (后面的 make_ref_struct, make_recon_struct, compute_nll_components 函数保持不变)
%% ======== 辅助函数 ========
function R = make_ref_struct(data)
    row0 = data(2, :); rowT = data(end, :);
    R.t0 = row0(1); R.mu0 = row0(3); 
    % 使用 (P95-P5)/4 作为 sigma 的估计
    R.sig0 = (row0(4) - row0(2)) / 4; 
    
    R.tT = rowT(1); R.muT = rowT(3);
    R.sigT = (rowT(4) - rowT(2)) / 4;
end

function R = make_recon_struct(T, i0, iT, meanName, p5Name, p95Name)
    R.mu0 = T.(meanName)(i0);
    R.sig0 = (T.(p95Name)(i0) - T.(p5Name)(i0)) / 4;
    
    R.muT = T.(meanName)(iT);
    R.sigT = (T.(p95Name)(iT) - T.(p5Name)(iT)) / 4;
end

%% ======== 核心：计算 NLL 分量 (Delta 和 Mu 分开) ========
function [nll_delta, nll_mu] = compute_nll_components(ref, rec)
    % 1. Delta (变化量) 部分
    delta_ref = ref.muT - ref.mu0;
    delta_rec = rec.muT - rec.mu0;
    
    % 误差传播: sigma_delta^2 = sigma_T^2 + sigma_0^2
    var_delta_total = (ref.sigT^2 + ref.sig0^2) + (rec.sigT^2 + rec.sig0^2);
    
    diff_delta = delta_rec - delta_ref;
    
    % NLL_Delta = 0.5 * [ ln(var) + diff^2 / var ]
    % 第一项惩罚不确定性过大，第二项惩罚偏离过大
    nll_delta = 0.5 * (log(var_delta_total) + (diff_delta^2 / var_delta_total));
    
    % 2. Mu (均值) 部分
    % 这里我们比较起始点 mu0 和 结束点 muT 的平均偏差
    % 或者更简单，比较整个过程的“中心点”位置偏差
    mu_ref_center = 0.5 * (ref.mu0 + ref.muT);
    mu_rec_center = 0.5 * (rec.mu0 + rec.muT);
    
    % 误差传播: sigma_mu^2 = 0.25 * (sigma_T^2 + sigma_0^2)
    var_mu_total = 0.25 * ((ref.sigT^2 + ref.sig0^2) + (rec.sigT^2 + rec.sig0^2));
    
    diff_mu = mu_rec_center - mu_ref_center;
    
    nll_mu = 0.5 * (log(var_mu_total) + (diff_mu^2 / var_mu_total));
    
    % 处理可能的 NaN
    if isnan(nll_delta), nll_delta = Inf; end
    if isnan(nll_mu), nll_mu = Inf; end
end