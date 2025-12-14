%%%%%%  5 box ocean model Ruoyuan Qiu 2025 modified from Zhao et al.,2023

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Define parameters   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CHEES_circulation(d13C_CO2_val, d13C_methane_val)
% 5-box ocean model 主函数
%
% 用法（在 MATLAB 里）：
%   CHEES_main(-5, -37)
%
% 如果命令行不传参数，就用默认值：
if nargin < 1
    d13C_CO2_val = -5;
end
if nargin < 2
    d13C_methane_val = -37;
end
clc;close all
outdir = sprintf('0.8_%g_%g', d13C_CO2_val, d13C_methane_val);
if ~exist(outdir, 'dir')
    mkdir(outdir);
end
%%%% === load data ===
d13C_p  = load('data/d13Cp_data.txt');  % [time, value]
SST     = load('data/SST.txt');        
pco2    = load('data/CO2.txt');          

%%%% === 参数扫描设置 ===
mu_vals      = -max(abs(d13C_CO2_val), 22) : -1 : -min(abs(d13C_methane_val), 35);      % [-22, -23, ..., -35] 共 14 个
TotalC_vals  = 4200:100:6000;       % [4200, 4300, ..., 6000] 共 19 个
TotalC_ref   = 5006;                % 5000 对应原来 52e12/139e12 的标度
sigma        = 4/1.644853;          % 截断正态的 sigma 保持不变

for imu = 1:numel(mu_vals)
    for iF = 1:numel(TotalC_vals)

        mu     = mu_vals(imu);
        TotalC = TotalC_vals(iF);
        scaleC = TotalC / TotalC_ref;

        fprintf('\n=== Running ensemble for mu = %.0f, TotalC = %d ===\n', mu, TotalC);

        num_runs = 50 ;
        results  = cell(num_runs, 9);

        for ii = 1: num_runs
            rng(ii);
            try
                disp(['  Ensemble member : ',num2str(ii)]);

                %%%%%%% set up global structures
                global stepnumber
                global pars
                global workingstate

                % --------- 这里是你原来的参数设置，只改了 d13C_input 和 input_C ---------
                %pars.d13C_input = [-25, mu, mu, mu];
                %pars.time       = [-3e6 0 100  6000];
                %pars.input_C    = [ 0  0 52e12 139e12] * scaleC;

                pars.time = [ ...
                    -3.00e+06 0.00e+00 0.1 5.00e+01 1.50e+02 2.50e+02 3.50e+02 4.50e+02 5.50e+02 ...
                    6.50e+02 7.50e+02 8.50e+02 9.50e+02 1.05e+03 1.15e+03 1.25e+03 1.35e+03 ...
                    1.45e+03 1.55e+03 1.65e+03 1.75e+03 1.85e+03 1.95e+03 2.05e+03 2.15e+03 ...
                    2.25e+03 2.35e+03 2.45e+03 2.55e+03 2.65e+03 2.75e+03 2.85e+03 2.95e+03 ...
                    3.05e+03 3.15e+03 3.25e+03 3.35e+03 3.45e+03 3.55e+03 3.65e+03 3.75e+03 ...
                    3.85e+03 3.95e+03 4.05e+03 4.15e+03 4.25e+03 4.35e+03 4.45e+03 ...
                    4.55e+03 4.65e+03 4.75e+03 4.85e+03 4.95e+03 5.05e+03 5.15e+03 5.25e+03 ...
                    5.35e+03 5.45e+03 5.55e+03 5.65e+03 5.75e+03 5.85e+03 5.95e+03 6000];

                % 2. 定义碳排放输入 (spin-up 阶段为 0，后续为图片中的数据，最后乘以 scaleC)
                pars.input_C = [ ...
                    0.00e+00 0.00e+00 2.52e13 2.72e+13 3.12e+13 3.42e+13 3.66e+13 3.81e+13 4.00e+13 ...
                    4.16e+13 4.31e+13 4.45e+13 4.55e+13 4.69e+13 4.77e+13 4.88e+13 4.99e+13 ...
                    5.06e+13 5.16e+13 5.25e+13 5.36e+13 5.46e+13 5.57e+13 5.66e+13 5.78e+13 ...
                    5.87e+13 5.97e+13 6.06e+13 6.19e+13 6.31e+13 6.44e+13 6.55e+13 6.67e+13 ...
                    6.78e+13 6.93e+13 7.03e+13 7.14e+13 7.29e+13 7.39e+13 7.48e+13 7.58e+13 ...
                    7.82e+13 7.97e+13 8.11e+13 8.24e+13 8.37e+13 8.52e+13 8.68e+13 8.83e+13 ...
                    8.98e+13 9.13e+13 9.31e+13 9.43e+13 9.65e+13 9.83e+13 1.00e+14 1.02e+14 ...
                    1.04e+14 1.06e+14 1.08e+14 1.10e+14 1.12e+14 1.13e+14 1.135e14] * scaleC;

                % 3. 定义同位素输入 (spin-up 开始时为 -25，一旦进入 t=0 及事件发生期全部设为 mu)
                % 长度必须与 pars.time 一致 (共 62 个点：2个 spin-up + 60个数据)
                pars.d13C_input = [-25, ones(1, 63) * mu];

                pars.k_ccdeg = 32.8e12 * (0.75+0.5*rand(1)) ; 
                pars.k_phosw = 2.5e+11 * (0.75+0.5*rand(1)) ;  

                %pars.d13C_CO2     = -5;      % 固定不变
                %pars.d13C_methane = -37;     % 固定不变
                pars.d13C_CO2     = d13C_CO2_val;      % 固定不变
                pars.d13C_methane = d13C_methane_val;     % 固定不变
                pars.methane_oxygen = 8.5e3 + 2.5e3 * rand(1);

                pars.k_CH4_diff          = 0.7 + 0.2 * rand(1);          
                pars.k_CH4_from_ocean    = 0.8  + 0.2 * rand(1);   
                pars.k_CH4_oxiwatercolum = 0.2 + 0.1 * rand(1); 

                %%%%%% water reservoir sizes in m3 (m=margins, s=surface, h= hi-lat, d=deep)
                pars.vol_p  = 2.6e15 ;  
                pars.vol_di = 5.4e15 ;  
                pars.vol_s  = 2.75e16 ; 
                pars.vol_h  = 1.22e16 ; 
                pars.vol_d  = 1.35e18 ;
                pars.vol_ocean = pars.vol_p + pars.vol_di + pars.vol_s + pars.vol_h + pars.vol_d ;

                %%%%%% mixing coefficient (Sv)
                pars.mixcoeff_dip = 30.28;
                pars.mixcoeff_ds  = 46.33;
                pars.mixcoeff_dh  = 54.9;

                %%%%%% inorganic carbon reservoirs in moles C
                pars.CO2_a_0  = 5e16 ; 
                pars.DIC_p_0  = 5.2e15 ;
                pars.DIC_di_0 = 1.08e16 ;
                pars.DIC_s_0  = 5.37e16 ;
                pars.DIC_h_0  = 2.71e16 ;
                pars.DIC_d_0  = 3e18 ;
                pars.ALK_p_0  = 5.2e15 ;
                pars.ALK_di_0 = 1.08e16 ;
                pars.ALK_s_0  = 5.37e16 ;
                pars.ALK_h_0  = 2.71e16 ;
                pars.ALK_d_0  = 3e18 ;

                %%%%%% C isotope composition
                pars.d13c_atm_0    = -7 ;
                pars.d13c_DIC_p_0  = 0.1 ;
                pars.d13c_DIC_di_0 = 0.1 ;
                pars.d13c_DIC_s_0  = 0.1 ;
                pars.d13c_DIC_h_0  = 0.1 ;
                pars.d13c_DIC_d_0  = 0.1 ;

                %%%%%% POC reservoirs in moles C
                pars.POC_p_0  = 630e12 ;
                pars.POC_di_0 = 250e12 ;
                pars.POC_s_0  = 2329e12 ;
                pars.POC_h_0  = 1084e12 ;
                pars.POC_d_0  = 56000e12 ;

                %%%%%% Dissolved phosphate in moles
                pars.DP_p_0  = 1.82e12 ;
                pars.DP_di_0 = 7.56e12 ;
                pars.DP_s_0  = 0.55e12 ;
                pars.DP_h_0  = 16.27e12 ;
                pars.DP_d_0  = 2970e12 ;

                pars.d13c_POC_p_0  = -26; 
                pars.d13c_POC_di_0 = -26;
                pars.d13c_POC_s_0  = -26 ;
                pars.d13c_POC_h_0  = -26 ;
                pars.d13c_POC_d_0  = -26 ;

                %%%%%% initial amount of O2 in moles C
                pars.O2_a_0  = 3.7e19 ; 
                pars.O2_p_0  = 6.705e14 ;
                pars.O2_di_0 = 8.964e14 ;
                pars.O2_s_0  = 9.139e15 ;
                pars.O2_h_0  = 4.02e15 ;
                pars.O2_d_0  = 1.823e17 ;

                pars.O2_conc_p_0  = pars.O2_p_0/pars.vol_p ;
                pars.O2_conc_di_0 = pars.O2_di_0/pars.vol_di;
                pars.O2_conc_s_0  = pars.O2_s_0/pars.vol_s;
                pars.O2_conc_h_0  = pars.O2_h_0/pars.vol_h;
                pars.O2_conc_d_0  = pars.O2_d_0/pars.vol_d;
                pars.Sredu_part   = 0.0860;

                %%%%%% initial amount of FeIII in moles
                pars.FeIII_p_0  = 1.56e9 ;
                pars.FeIII_di_0 = 3.24e9 ;
                pars.FeIII_s_0  = 9.625e9 ;
                pars.FeIII_h_0  = 3.66e9 ;
                pars.FeIII_d_0  = 810e9 ;

                %%%%%% initial amount of sulfate in moles
                pars.SO4_p_0  = 3.12e16;
                pars.SO4_di_0 = 5.4e15*12;
                pars.SO4_s_0  = 2.75e16*12;
                pars.SO4_h_0  = 1.22e16 *12;
                pars.SO4_d_0  = 1.35e18 *12;

                %%%%%% initial amount of FeII in moles
                pars.FeII_p_0  = 0;
                pars.FeII_di_0 = 0;
                pars.FeII_s_0  = 0;
                pars.FeII_h_0  = 0;
                pars.FeII_d_0  = 0;

                %%%%%% initial amount of H2S in moles
                pars.H2S_p_0  = 0;
                pars.H2S_di_0 = 0;
                pars.H2S_s_0  = 0;
                pars.H2S_h_0  = 0;
                pars.H2S_d_0  = 0;

                %%%%%% present day rates
                pars.k_carbw = 12e12 ; 
                pars.k_sfw   = 0 ;       
                pars.k_mccb  = 20e12;  
                pars.k_silw  = pars.k_mccb - pars.k_carbw ;  
                basfrac      = 0.3 ;

                pars.k_granw = pars.k_silw * (1 - basfrac) ;
                pars.k_basw  = pars.k_silw * basfrac ;

                %%%%%% organic C cycle
                pars.k_ocdeg = 1.25e12 * 2  ;
                pars.k_locb  = 2.5e12 ; 
                pars.k_mocb  = 7e12  ;  
                pars.k_oxidw = pars.k_mocb + pars.k_locb - pars.k_ocdeg ;

                %%%%%% present P, Fe, pyrite and sulfate weathering rate
                pars.k_FeIIIw   = 2e9;
                pars.k_pyritew  = 1.85e12;
                pars.k_sulfatew = 1.25e12;

                %%%%%% present pyrite and sulfate burial rate
                pars.k_pyriteb_p  = 1.4e12;
                pars.k_pyriteb_di = 0.45e12;
                pars.k_sulfateb   = 1.25e12;

                %%%%%% Redfeild ratio
                pars.Red_C_P = 106;
                pars.Red_C_N = 106/16;
                pars.Red_C_O = 106/138;
                pars.Red_C_Fe = 106*2000;
                pars.Red_Fe_P = pars.Red_C_P/pars.Red_C_Fe;

                pars.BE_p  = 6/40;
                pars.BE_di = 6/40;
                pars.BE_d  = 1/10;

                %%%%%% Monod constant; mol/m3
                pars.KP      = 0.1e-3;
                pars.KFe     = 0.1e-6;
                pars.KmO2    = 10e-3;  
                pars.KmO2s   = 10e-3;  
                pars.KmO2m   = 1e-3;   
                pars.KmFeIII = 10;
                pars.KmSO4   = 0.5;    

                %%%%%% reaction rate constant m3 mol-1 yr-1
                pars.kpy     = 0.3708/1e3*24*365.25;
                pars.kironO  = 1.4e5;
                pars.ksulfO  = 1.6e2; 
                pars.kSironR = 8;
                pars.kSide   = 4e3;   

                %%%%%% Ksp
                pars.Kspside   = 10^(-8.4)*1e6;
                pars.KspFeSaq  = 10^(-5.08)*1e3;
                pars.STFeSaq   = 10^(-5.7)*1e3;

                pars.FeIIIa_s  = 1.443e9; 
                pars.FeIIIa_h  = 0;  
                pars.FeIIhydro = 13.5e9;

                %%%%%% initial amount of CH4 in moles C
                pars.CH4_a_0      = 8e13; 
                pars.f_CH4_a      = 10e12; 
                pars.methaneHydrate = 5.25e17; 

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%   Initialise   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                options = odeset('maxstep',1e3,'NonNegative',[1:11 18:27 33:59]) ;

                stepnumber = 1 ;

                pars.whenstart = -3e6 ;
                pars.whenend   = 6000 ; 

                %%%%% model start state
                pars.startstate(1)  = pars.CO2_a_0 ;
                pars.startstate(2)  = pars.DIC_p_0 ;
                pars.startstate(3)  = pars.DIC_di_0 ;
                pars.startstate(4)  = pars.DIC_s_0 ;
                pars.startstate(5)  = pars.DIC_h_0 ;
                pars.startstate(6)  = pars.DIC_d_0 ;

                pars.startstate(7)  = pars.ALK_p_0 ;
                pars.startstate(8)  = pars.ALK_di_0 ;
                pars.startstate(9)  = pars.ALK_s_0 ;
                pars.startstate(10) = pars.ALK_h_0 ;
                pars.startstate(11) = pars.ALK_d_0 ;

                pars.startstate(12) = pars.CO2_a_0 * pars.d13c_atm_0 ;
                pars.startstate(13) = pars.DIC_p_0  * pars.d13c_DIC_p_0 ;
                pars.startstate(14) = pars.DIC_di_0 * pars.d13c_DIC_di_0 ;
                pars.startstate(15) = pars.DIC_s_0  * pars.d13c_DIC_s_0 ;
                pars.startstate(16) = pars.DIC_h_0  * pars.d13c_DIC_h_0 ;
                pars.startstate(17) = pars.DIC_d_0  * pars.d13c_DIC_d_0 ;

                pars.startstate(18) = pars.POC_p_0;
                pars.startstate(19) = pars.POC_di_0;
                pars.startstate(20) = pars.POC_s_0;
                pars.startstate(21) = pars.POC_h_0;
                pars.startstate(22) = pars.POC_d_0;

                pars.startstate(23) = pars.DP_p_0;
                pars.startstate(24) = pars.DP_di_0;
                pars.startstate(25) = pars.DP_s_0;
                pars.startstate(26) = pars.DP_h_0;
                pars.startstate(27) = pars.DP_d_0;

                pars.startstate(28) = pars.d13c_POC_p_0  *pars.POC_p_0;
                pars.startstate(29) = pars.d13c_POC_di_0 *pars.POC_di_0;
                pars.startstate(30) = pars.d13c_POC_s_0  *pars.POC_s_0;
                pars.startstate(31) = pars.d13c_POC_h_0  *pars.POC_h_0;
                pars.startstate(32) = pars.d13c_POC_d_0  *pars.POC_d_0;

                pars.startstate(33) = pars.O2_a_0;
                pars.startstate(34) = pars.O2_p_0;
                pars.startstate(35) = pars.O2_di_0;
                pars.startstate(36) = pars.O2_s_0;
                pars.startstate(37) = pars.O2_h_0;
                pars.startstate(38) = pars.O2_d_0;

                pars.startstate(39) = pars.FeIII_p_0;
                pars.startstate(40) = pars.FeIII_di_0;
                pars.startstate(41) = pars.FeIII_s_0;
                pars.startstate(42) = pars.FeIII_h_0;
                pars.startstate(43) = pars.FeIII_d_0;

                pars.startstate(44) = pars.SO4_p_0;
                pars.startstate(45) = pars.SO4_di_0;
                pars.startstate(46) = pars.SO4_s_0;
                pars.startstate(47) = pars.SO4_h_0;
                pars.startstate(48) = pars.SO4_d_0;

                pars.startstate(49) = pars.FeII_p_0;
                pars.startstate(50) = pars.FeII_di_0;
                pars.startstate(51) = pars.FeII_s_0;
                pars.startstate(52) = pars.FeII_h_0;
                pars.startstate(53) = pars.FeII_d_0;

                pars.startstate(54) = pars.H2S_p_0;
                pars.startstate(55) = pars.H2S_di_0;
                pars.startstate(56) = pars.H2S_s_0;
                pars.startstate(57) = pars.H2S_h_0;
                pars.startstate(58) = pars.H2S_d_0;

                pars.startstate(59) = pars.CH4_a_0;
                pars.startstate(60) = pars.methaneHydrate;

                %%%%%%% run the system  
                [rawoutput.T,rawoutput.Y] = ode15s(@CHEES_equations,[pars.whenstart pars.whenend],pars.startstate,options);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%   Postprocessing   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                pars.output_length = length(rawoutput.T) ;

                [sharedvals,trecords] = intersect(workingstate.time,rawoutput.T,'stable') ;

                field_names = fieldnames(workingstate) ;
                for numfields = 1:length(field_names)
                    eval([' state.' char( field_names(numfields) ) ' = workingstate.' char( field_names(numfields) ) '(trecords) ; '])
                end

                results{ii, 1} = pars.k_ccdeg; 
                results{ii, 2} = pars.k_phosw;          
                results{ii, 3} = pars.methane_oxygen;
                results{ii, 4} = pars.k_CH4_diff;
                results{ii, 5} = pars.k_CH4_from_ocean;
                results{ii, 6} = pars.k_CH4_oxiwatercolum;
                results{ii, 7} = pars.k_CH4_oxiwatercolum;
                results{ii, 8} = pars.k_CH4_oxiwatercolum;
                results{ii, 9} = state;     

            catch ME
                fprintf('Error in run %d: %s\n', ii, ME.message);
                results{ii, 1} = [];
                results{ii, 2} = [];
                results{ii, 3} = [];
                results{ii, 4} = [];
                results{ii, 5} = [];
                results{ii, 6} = [];
                results{ii, 7} = [];
                results{ii, 8} = [];
                results{ii, 9} = [];
            end

        end % for ii = 1:num_runs

        %%%%%%%%%%%%%% ensemble 后处理 %%%%%%%%%%%%%%
        % ... (for ii = 1:num_runs 循环结束) ...

        %% %%%%%%%%%%%%%% ensemble 后处理 (Robust Version) %%%%%%%%%%%%%%
        
        % =========================================================================
    % ================== 以下为后处理核心修正部分 (请完全替换) ==================
    % =========================================================================

    %% 1. 剔除失败的 Run (防止空数据导致后续报错)
    % 移除完全为空的行
    empty_rows = all(cellfun(@(x) isempty(x) || (isstring(x) && strlength(x)==0), results), 2);
    results = results(~empty_rows, :);

    % 【关键】进一步检查：移除没有运行到 t >= 0 的 Run
    % 如果 ode15s 在负年份就崩了(也就是你日志里的警告)，这些数据会导致后续插值失败
    valid_mask = false(size(results, 1), 1);
    for k = 1:size(results, 1)
        try
            if isstruct(results{k, 9}) && isfield(results{k, 9}, 'time_myr') 
                % 只要最大时间能覆盖到 0 以后，就算成功
                if max(results{k, 9}.time_myr) >= 0
                    valid_mask(k) = true;
                end
            end
        catch
            valid_mask(k) = false;
        end
    end
    results = results(valid_mask, :);
    
    num_runs_eff = size(results, 1);
    
    % 如果这一组参数导致所有 Run 都挂了，直接跳过，不要报错退出
    if num_runs_eff == 0
        fprintf('  [Warning] All runs failed for mu=%.1f, TotalC=%.1f. Skipping.\n', mu, TotalC);
        continue; % 跳过当前 TotalC 循环，进入下一组
    end

    %% 2. 建立统一时间轴并插值
    
    % 找到最长的一条时间轴作为基准
    time_lengths = zeros(num_runs_eff, 1);
    for i = 1:num_runs_eff
        time_lengths(i) = length(results{i, 9}.time_myr); 
    end
    [~, max_idx] = max(time_lengths); 
    base_time = results{max_idx, 9}.time_myr;
    
    % 强制 base_time 为列向量 (N x 1)
    base_time = base_time(:);

    % 只有当最长的运行时间太短时（比如都没跑到6000年），才强制使用标准轴，防止插值范围不足
    if max(base_time) < 0.006 
         base_time = pars.time(pars.time >= -100 & pars.time <= 6000)' / 1e6; 
         base_time = base_time(:);
    end

    len_t = length(base_time);
    
    % 初始化大矩阵 (行=Run, 列=Time)
    d13c_DIC            = NaN(num_runs_eff, len_t);
    T_p                 = NaN(num_runs_eff, len_t);
    Atmospheric_CO2_ppm = NaN(num_runs_eff, len_t);
    Atmospheric_CH4_a_ppm = NaN(num_runs_eff, len_t);
    f_volcanic_CO2      = NaN(num_runs_eff, len_t);
    f_CH4               = NaN(num_runs_eff, len_t);
    d13C_input          = NaN(num_runs_eff, len_t);

    for i = 1:num_runs_eff
        try
            state_data   = results{i, 9};
            time_myr_run = state_data.time_myr;
            
            % 定义插值函数: linear插值, 超出范围填NaN
            do_interp = @(vec) interp1(time_myr_run, vec, base_time, 'linear', NaN);

            % === 提取 DIC 相关并计算加权 ===
            M_p = do_interp(state_data.DIC_p);
            M_s = do_interp(state_data.DIC_s);
            M_h = do_interp(state_data.DIC_h);
            D_p = do_interp(state_data.d13c_DIC_p);
            D_s = do_interp(state_data.d13c_DIC_s);
            D_h = do_interp(state_data.d13c_DIC_h);
            
            total_mass = M_p + M_s + M_h;
            weighted_d13c = (M_p .* D_p + M_s .* D_s + M_h .* D_h) ./ total_mass;
            
            % 【关键修正】赋值前强制转为行向量 (1 x N)，防止维度不匹配报错
            d13c_DIC(i, :) = weighted_d13c(:)'; 
            
            % === 提取其他变量 (同样强制转行向量) ===
            vec_Tp   = do_interp(state_data.T_p);
            T_p(i, :) = vec_Tp(:)';
            
            vec_CO2  = do_interp(state_data.Atmospheric_CO2_ppm);
            Atmospheric_CO2_ppm(i, :) = vec_CO2(:)';
            
            vec_CH4a = do_interp(state_data.Atmospheric_CH4_a_ppm);
            Atmospheric_CH4_a_ppm(i,:) = vec_CH4a(:)';
            
            vec_fVolc = do_interp(state_data.f_volcanic_CO2);
            f_volcanic_CO2(i, :) = vec_fVolc(:)';
            
            vec_fCH4 = do_interp(state_data.f_CH4);
            f_CH4(i, :) = vec_fCH4(:)';
            
            % 注意：d13c_CO2_input 可能在某些版本代码中叫 pars.d13C_input
            % 如果 state 里没有这个变量，这里会 catch 到错误并保持 NaN，不会崩
            if isfield(state_data, 'd13c_CO2_input')
                vec_d13Cin = do_interp(state_data.d13c_CO2_input);
                d13C_input(i, :) = vec_d13Cin(:)';
            end
            
        catch ME
            % 仅仅打印警告，不要中断整个循环
            % fprintf('  [Warning] Interpolation failed for run index %d\n', i);
        end
    end

    %% 3. 统计计算 (忽略 NaN)
    
    % 定义辅助函数：同时算 Mean, 5%, 95%
    % 使用 'omitnan' 确保即使有部分 NaN 也能算出均值
    calc_stats = @(mat) deal(mean(mat, 1, 'omitnan'), prctile(mat, 5, 1), prctile(mat, 95, 1));

    [d13C_input_mean, d13C_input_5, d13C_input_95]    = calc_stats(d13C_input);
    [f_volcanic_CO2_mean, f_volcanic_CO2_5, f_volcanic_CO2_95] = calc_stats(f_volcanic_CO2);
    [f_CH4_mean, f_CH4_5, f_CH4_95]                   = calc_stats(f_CH4);
    [d13c_DIC_mean, d13c_DIC_5, d13c_DIC_95]          = calc_stats(d13c_DIC);
    [T_p_mean, T_p_5, T_p_95]                         = calc_stats(T_p);
    [Atmospheric_CO2_mean, Atmospheric_CO2_5, Atmospheric_CO2_95] = calc_stats(Atmospheric_CO2_ppm);
    [Atmospheric_CH4_mean, Atmospheric_CH4_5, Atmospheric_CH4_95] = calc_stats(Atmospheric_CH4_a_ppm);

    %% 4. 生成 CSV (修复 table 维度报错)
    
    plot_time = base_time * 1e6; 
    
    % 只保留 t >= 0 的时间点
    idx = plot_time >= 0; 
    
    % 如果没有正半轴数据，跳过写入
    if sum(idx) == 0
         continue;
    end

    % 【终极修复】强制列向量转换函数
    % 无论输入是行还是列，reshape(..., [], 1) 都会把它变成列向量
    % 这就恢复了你原来代码中 val = val(:) 的那种稳健性
    make_col = @(vec) reshape(vec(idx), [], 1); 

    % 时间列
    time_yr_col = make_col(plot_time);

    % 组装 Table
    % 这里的每一个输入都经过 make_col 处理，保证行数绝对一致
    TS_out = table( ...
        time_yr_col, ...
        make_col(d13c_DIC_mean), make_col(d13c_DIC_5), make_col(d13c_DIC_95), ...
        make_col(Atmospheric_CO2_mean), make_col(Atmospheric_CO2_5), make_col(Atmospheric_CO2_95), ...
        make_col(T_p_mean - 273), make_col(T_p_5 - 273), make_col(T_p_95 - 273), ... % 转为摄氏度
        make_col(f_CH4_mean), make_col(f_CH4_5), make_col(f_CH4_95), ...
        make_col(f_volcanic_CO2_mean), make_col(f_volcanic_CO2_5), make_col(f_volcanic_CO2_95), ...
        make_col(d13C_input_mean), make_col(d13C_input_5), make_col(d13C_input_95), ...
        make_col(Atmospheric_CH4_mean), make_col(Atmospheric_CH4_5), make_col(Atmospheric_CH4_95), ...
        'VariableNames', { ...
            'time_yr', ...
            'd13C_DIC_mean', 'd13C_DIC_p5', 'd13C_DIC_p95', ...
            'CO2_mean_ppm',  'CO2_p5_ppm',  'CO2_p95_ppm', ...
            'SST_mean_degC', 'SST_p5_degC', 'SST_p95_degC', ...
            'f_CH4_mean_mol_per_yr',        'f_CH4_p5_mol_per_yr',        'f_CH4_p95_mol_per_yr', ...
            'f_volcanic_CO2_mean_mol_per_yr','f_volcCO2_p5_mol_per_yr',   'f_volcCO2_p95_mol_per_yr', ...
            'd13C_input_mean', 'd13C_input_p5', 'd13C_input_p95', ...
            'CH4_mean_ppm', 'CH4_p5_ppm', 'CH4_p95_ppm' ...
        } ...
    );

    filename = sprintf('CHEES_%d_%d.csv', round(mu), round(TotalC));
    fullpath = fullfile(outdir, filename);

    try
        writetable(TS_out, fullpath);
        fprintf('  Saved %s\n', fullpath);
    catch ME
        fprintf('  [Error] Failed to write CSV: %s\n', ME.message);
    end


    end % for iF
end % for imu
end

%% ---- Local functions (must be at end of script) ----
function x = truncnorm(mu, sigma, lo, hi)
x = mu + sigma*randn();
while x < lo || x > hi
    x = mu + sigma*randn();
end
end
