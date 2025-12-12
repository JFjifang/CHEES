% 
% function [total_ll,ll_d13C, ll_CO2, ll_SST] = log_likelihood(time_myr,d13C,d13c_POC,pco2,Atmospheric_CO2,SST,GAST)
% 
% 
%     [~, idx0] = min(abs(time_myr*1e-6 - 0));   % 或 min(abs(time))
%     time_at0  = time_myr(idx0);    
%     d13C_at0  = d13c_POC(idx0);
%     Atmospheric_CO2_at0  = Atmospheric_CO2(idx0);
%     GAST_at0  = GAST(idx0)-273;
% 
%     [~, idx6000] = min(abs(time_myr*1e-6 - 6000));
%     time_at6000  = time_myr(idx6000);          
%     d13C_at6000  = d13c_POC(idx6000);
%     Atmospheric_CO2_at6000  = Atmospheric_CO2(idx6000);
%     GAST_at6000  = GAST(idx6000)-273;
% 
%     d13C_model = [d13C_at0;d13C_at6000];
%     CO2_model  = [Atmospheric_CO2_at0;Atmospheric_CO2_at6000];
%     SST_model  = [GAST_at0;GAST_at6000];
% 
%     ll_d13C = sum(log(pdf('Normal', d13C(:,3), d13C_model, 0.15)));  % d13C(1,3) - d13C_at0 +    d13C(2,3) - d13C_at6000 ; 
%     ll_CO2  = sum(log(pdf('Normal', pco2(:,3)./100, CO2_model./100, 0.15))); % (pco2(1,3) - Atmospheric_CO2_at0 +    pco2(2,3) - Atmospheric_CO2_at6000) / 1000 ; 
%     ll_SST  = sum(log(pdf('Normal', SST(:,3)./10, SST_model./10, 0.15))); %(SST(1,3) - (GAST_at0-273) +    SST(2,3) - (GAST_at6000-273)) / 10 ; 
% 
%     total_ll = ll_d13C  + ll_CO2 + ll_SST;
% 
% 
% 
% end
% 
function [total_ll, ll_d13C, ll_CO2, ll_SST] = log_likelihood(time_myr, d13C_data, d13c_POC_model, pco2_data, Atmospheric_CO2_model, SST_data, GAST_model)

    % ================= 配置区域 =================
    k_weight = 3.0; % 权重：依然保留，用于平衡“趋势”和“绝对值”
    
    % 注意：不再需要人为设定 sigma_d13C, sigma_CO2 等
    
    % ================= 数据提取与 Sigma 计算 =================
    
    % --- 1. d13C 数据 ---
    % 第1行对应 t=0 (Start)，最后一行对应 t=6000 (End)
    obs_d13C_0 = d13C_data(1, 3);
    obs_d13C_T = d13C_data(end, 3);
    
    % 动态计算 Sigma: (High - Low) / 4
    sigma_d13C_0 = (d13C_data(1, 4) - d13C_data(1, 2)) / 4;
    sigma_d13C_T = (d13C_data(end, 4) - d13C_data(end, 2)) / 4;
    
    % 安全检查：防止 sigma 为 0 (如果 High=Low)
    if sigma_d13C_0 < 1e-6, sigma_d13C_0 = 0.15; end 
    if sigma_d13C_T < 1e-6, sigma_d13C_T = 0.15; end


    % --- 2. CO2 数据 ---
    obs_CO2_0 = pco2_data(1, 3); % 直接读取 ppm
    obs_CO2_T = pco2_data(end, 3);
    
    sigma_CO2_0 = (pco2_data(1, 4) - pco2_data(1, 2)) / 4;
    sigma_CO2_T = (pco2_data(end, 4) - pco2_data(end, 2)) / 4;
    
    if sigma_CO2_0 < 1e-6, sigma_CO2_0 = 30.0; end % 防止除以0，给个保底值
    if sigma_CO2_T < 1e-6, sigma_CO2_T = 30.0; end


    % --- 3. SST 数据 ---
    obs_SST_0 = SST_data(1, 3);
    obs_SST_T = SST_data(end, 3);
    
    sigma_SST_0 = (SST_data(1, 4) - SST_data(1, 2)) / 4;
    sigma_SST_T = (SST_data(end, 4) - SST_data(end, 2)) / 4;
    
    if sigma_SST_0 < 1e-6, sigma_SST_0 = 2.0; end
    if sigma_SST_T < 1e-6, sigma_SST_T = 2.0; end


    % ================= 模型数据提取 =================
    [~, idx0] = min(abs(time_myr*1e-6 - 0));
    [~, idxT] = min(abs(time_myr*1e-6 - 6000));
    
    mod_d13C_0 = d13c_POC_model(idx0);
    mod_d13C_T = d13c_POC_model(idxT);
    
    mod_CO2_0  = Atmospheric_CO2_model(idx0);
    mod_CO2_T  = Atmospheric_CO2_model(idxT);
    
    mod_SST_0  = GAST_model(idx0) - 273; 
    mod_SST_T  = GAST_model(idxT) - 273;


    % ================= 计算 NLL =================
    % 注意：现在我们需要传入两个 sigma (sigma_0 和 sigma_T)
    
    [nll_d13C_delta, nll_d13C_mu] = calculate_component_nll(obs_d13C_0, obs_d13C_T, mod_d13C_0, mod_d13C_T, sigma_d13C_0, sigma_d13C_T);
    
    [nll_CO2_delta, nll_CO2_mu]   = calculate_component_nll(obs_CO2_0, obs_CO2_T, mod_CO2_0, mod_CO2_T, sigma_CO2_0, sigma_CO2_T);
    
    [nll_SST_delta, nll_SST_mu]   = calculate_component_nll(obs_SST_0, obs_SST_T, mod_SST_0, mod_SST_T, sigma_SST_0, sigma_SST_T);


    % ================= 加权合并 =================
    score_d13C = nll_d13C_delta + (1 / k_weight^2) * nll_d13C_mu;
    score_CO2  = nll_CO2_delta  + (1 / k_weight^2) * nll_CO2_mu;
    score_SST  = nll_SST_delta  + (1 / k_weight^2) * nll_SST_mu;
    
    ll_d13C = -score_d13C;
    ll_CO2  = -score_CO2;
    ll_SST  = -score_SST;
    
    total_ll = ll_d13C + ll_CO2 + ll_SST;
end

% 辅助函数：更新为支持不同的 sigma_0 和 sigma_T
function [nll_delta, nll_mu] = calculate_component_nll(ref_0, ref_T, rec_0, rec_T, sigma_ref_0, sigma_ref_T)
    sigma_rec = 0; % 模型无误差
    
    % --- Delta (趋势) ---
    delta_ref = ref_T - ref_0;
    delta_rec = rec_T - rec_0;
    
    % 误差传播：Var(Delta) = Var(T) + Var(0)
    var_delta = (sigma_ref_T^2 + sigma_ref_0^2) + (sigma_rec^2 + sigma_rec^2);
    
    nll_delta = 0.5 * (log(var_delta) + (delta_rec - delta_ref)^2 / var_delta);
    
    % --- Mu (均值) ---
    mu_ref = 0.5 * (ref_T + ref_0);
    mu_rec = 0.5 * (rec_T + rec_0);
    
    % 误差传播：Var(Mean) = 0.25 * (Var(T) + Var(0))
    var_mu = 0.25 * ((sigma_ref_T^2 + sigma_ref_0^2) + (sigma_rec^2 + sigma_rec^2));
    
    nll_mu = 0.5 * (log(var_mu) + (mu_rec - mu_ref)^2 / var_mu);
end