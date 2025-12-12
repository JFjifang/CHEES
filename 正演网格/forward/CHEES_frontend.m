
%%%%%%  5 box ocean model Ruoyuan Qiu 2025 modified from Zhao et al.,2023

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Define parameters   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CHEES_frontend(d13C_CO2_val, d13C_methane_val, mu_val, TotalC)
clc;close all
outdir = sprintf('%d_%d', round(d13C_CO2_val), round(d13C_methane_val));
if ~exist(outdir, 'dir')
    mkdir(outdir);
end
%%%% === load data ===
d13C_p  = load('data/d13Cp_data.txt');  % [time, value]
SST     = load('data/SST.txt');        
pco2    = load('data/CO2.txt');          

num_runs = 50 ;
results = cell(num_runs, 9);

for ii = 1: num_runs
    rng(ii);
    try
    disp(['Number : ',num2str(ii)]);

%%%%%%% set up global structures
global stepnumber
global pars
global workingstate

mu    = mu_val;
sigma = 4/1.644853;
TotalC_ref   = 5000;
scaleC = TotalC / TotalC_ref;
pars.d13C_input = [-25, mu, mu, mu];
pars.time       = [-3e6 0 100  6000];
pars.input_C    = [ 0  0 52e12 139e12] * scaleC; %90/190
%pars.d13C_input = [-25 -22-13*rand(1) -22-13*rand(1) -22-13*rand(1)];%-23-4*rand(1) -23-4*rand(1) -23-4*rand(1)
pars.k_ccdeg = 21e12 * (0.75+0.5*rand(1)) ; % 8e12 18e12 * (0.75+0.5*rand(1))
pars.k_phosw = 2.5e+11 * (0.75+0.5*rand(1)) ;  % 1.25 1.2088e+11 * (0.75+0.5*rand(1))

%pars.d13C_CO2 = -6;
%pars.d13C_methane = -37;%-60
pars.d13C_CO2     = d13C_CO2_val;      % 固定不变
pars.d13C_methane = d13C_methane_val;     % 固定不变
pars.methane_oxygen = 8.5e3 + 2.5e3 * rand(1);

pars.k_CH4_diff = 0.8 + 0.15 * rand(1);          % 有多少甲烷从沉积物中扩散到海洋中，原本是0.8，0.15
pars.k_CH4_from_ocean = 0.8  + 0.2 * rand(1);   % 有多少甲烷是来自海洋的0.4+ 0.3 * rand(1)
pars.k_CH4_oxiwatercolum = 0.2 + 0.1 * rand(1); % 这些甲烷气泡在水柱中有多少被氧化，0.1 + 0.3 * rand(1)

%%%%%% water reservoir sizes in m3 (m=margins, s=surface, h= hi-lat, d=deep)
pars.vol_p  = 2.6e15 ;  %%%% approx volume of all shelves and slope to depth 100m, area pecentage 5%
pars.vol_di = 5.4e15 ;  %%%% approx volume of all shelves and slope in depth 100-1000m, area pecentage 5%
pars.vol_s  = 2.75e16 ; %%%% approx volume of suface water to depth 100m, area pecentage 76.5%
pars.vol_h  = 1.22e16 ; %%%% approx volume of hi-lat to depth 250m, area pecentage 13.5%
pars.vol_d  = 1.35e18 ;
pars.vol_ocean = pars.vol_p + pars.vol_di + pars.vol_s + pars.vol_h + pars.vol_d ;

%%%% mixing coefficient (Sv)
pars.mixcoeff_dip = 30.28;
pars.mixcoeff_ds  = 46.33;
pars.mixcoeff_dh  = 54.9;

%%%%%% inorganic carbon reservoirs in moles C
pars.CO2_a_0  = 5e16 ;% default 5e16 
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

pars.d13c_POC_p_0  = -26; % -26
pars.d13c_POC_di_0 = -26;
pars.d13c_POC_s_0  = -26 ;
pars.d13c_POC_h_0  = -26 ;
pars.d13c_POC_d_0  = -26 ;

%%%%%% initial amount of O2 in moles C
pars.O2_a_0  = 3.7e19 ; % 3.7e19
pars.O2_p_0  = 6.705e14 ;
pars.O2_di_0 = 8.964e14 ;
pars.O2_s_0  = 9.139e15 ;
pars.O2_h_0  = 4.02e15 ;
pars.O2_d_0  = 1.823e17 ;

pars.O2_conc_p_0 = pars.O2_p_0/pars.vol_p ;
pars.O2_conc_di_0 =pars.O2_di_0/pars.vol_di;
pars.O2_conc_s_0 =pars.O2_s_0/pars.vol_s;
pars.O2_conc_h_0 =pars.O2_h_0/pars.vol_h;
pars.O2_conc_d_0 =pars.O2_d_0/pars.vol_d;
pars.Sredu_part = 0.0860;

%%%%%% initial amount of FeIII in moles
pars.FeIII_p_0  = 1.56e9 ;
pars.FeIII_di_0 = 3.24e9 ;
pars.FeIII_s_0  = 9.625e9 ;
pars.FeIII_h_0  = 3.66e9 ;
pars.FeIII_d_0  = 810e9 ;

%%%%%% initial amount of sulfate in moles
%%%% 12 mM
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

pars.k_carbw = 12e12 ; % 12e12
pars.k_sfw   = 0 ;       %%% Seafloor weathering
pars.k_mccb  = 20e12;  % pars.k_carbw + pars.k_ccdeg - pars.k_sfw  20e12
pars.k_silw  = pars.k_mccb - pars.k_carbw ;  % silicate weathering
basfrac = 0.3 ;

pars.k_granw = pars.k_silw * (1 - basfrac) ;
pars.k_basw = pars.k_silw * basfrac ;

%%%%%% organic C cycle
pars.k_ocdeg = 1.25e12 * 2  ;
pars.k_locb = 2.5e12 ; % 4.5e12 ; 2.5e12 
pars.k_mocb = 7e12  ;  % 4.5e12 ;
pars.k_oxidw = pars.k_mocb + pars.k_locb - pars.k_ocdeg ;

%%%%%% present P, Fe, pyrite and sulfate weathering rate
pars.k_FeIIIw = 2e9;
pars.k_pyritew = 1.85e12;
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
pars.KP = 0.1e-3;
pars.KFe = 0.1e-6;
pars.KmO2 = 10e-3;  %%%for sulfate reduction in water column
pars.KmO2s = 10e-3;  %%%for sulfate reduction in sediments 10e-3
pars.KmO2m = 1e-3;  %%%for methane oxidation in water column
pars.KmFeIII = 10;
pars.KmSO4 = 0.5; %%%0.5,1.6


%%%%%% reaction rate constant m3 mol-1 yr-1
pars.kpy = 0.3708/1e3*24*365.25;
pars.kironO = 1.4e5;
pars.ksulfO = 1.6e2; %%%1.6e2
pars.kSironR = 8;
pars.kSide = 4e3; %%%mol m-3 yr-1

%%%%%% Ksp
pars.Kspside = 10^(-8.4)*1e6;%%%mol2 m-6
pars.KspFeSaq = 10^(-5.08)*1e3;%%%mol m-3

pars.STFeSaq = 10^(-5.7)*1e3;%%%mol m-3

pars.FeIIIa_s = 1.443e9; %%%
pars.FeIIIa_h = 0;  %%%% 0.00795e9
pars.FeIIhydro = 13.5e9;

%%%%%% initial amount of CH4 in moles C
pars.CH4_a_0 =  8e13; %   
pars.f_CH4_a = 10e12; %  

pars.methaneHydrate = 5.25e17; % 1e4 Gt CH4  = 1e19 g = 6.25e17 mol

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Initialise   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% set maximum step size for solver
% options = odeset('maxstep',1e3) ;
options = odeset('maxstep',1e3,'NonNegative',[1:11 18:27 33:59]) ;

%%%% run beginning
% fprintf('Beginning run: \n')

%%%% set stepnumber to 1
stepnumber = 1 ;

%%%%%%% model timeframe in years (0 = present day)
pars.whenstart = -3e6 ;
pars.whenend = 6000 ; % 0.33e6 ;

%%%% model start state
pars.startstate(1) = pars.CO2_a_0 ;
pars.startstate(2) = pars.DIC_p_0 ;
pars.startstate(3) = pars.DIC_di_0 ;
pars.startstate(4) = pars.DIC_s_0 ;
pars.startstate(5) = pars.DIC_h_0 ;
pars.startstate(6) = pars.DIC_d_0 ;

pars.startstate(7) = pars.ALK_p_0 ;
pars.startstate(8) = pars.ALK_di_0 ;
pars.startstate(9) = pars.ALK_s_0 ;
pars.startstate(10) = pars.ALK_h_0 ;
pars.startstate(11) = pars.ALK_d_0 ;

pars.startstate(12) = pars.CO2_a_0 * pars.d13c_atm_0 ;
pars.startstate(13) = pars.DIC_p_0 * pars.d13c_DIC_p_0 ;
pars.startstate(14) = pars.DIC_di_0 * pars.d13c_DIC_di_0 ;
pars.startstate(15) = pars.DIC_s_0 * pars.d13c_DIC_s_0 ;
pars.startstate(16) = pars.DIC_h_0 * pars.d13c_DIC_h_0 ;
pars.startstate(17) = pars.DIC_d_0 * pars.d13c_DIC_d_0 ;

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

pars.startstate(28) = pars.d13c_POC_p_0*pars.POC_p_0;
pars.startstate(29) = pars.d13c_POC_di_0*pars.POC_di_0;
pars.startstate(30) = pars.d13c_POC_s_0*pars.POC_s_0;
pars.startstate(31) = pars.d13c_POC_h_0*pars.POC_h_0;
pars.startstate(32) = pars.d13c_POC_d_0*pars.POC_d_0;

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

%%%%% start time counter
% tic



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Postprocessing   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% takes 'workingstate' from model and turns into 'state' %%%%%%%

%%%% size of output 
pars.output_length = length(rawoutput.T) ;
%%%%%%%%%% model finished output to screen
% fprintf('Integration finished \t') ; fprintf('Total steps: %d \t' , stepnumber ) ; fprintf('Output steps: %d \n' , pars.output_length ) 
% toc

%%%%%%%%% print final model states using final state for each timepoint
%%%%%%%%% during integration
% fprintf('assembling state vectors... \t')
% tic

%%%% trecords is index of shared values between ode15s output T vector and
%%%% model recorded workingstate t vector
[sharedvals,trecords] = intersect(workingstate.time,rawoutput.T,'stable') ;

%%%%%% assemble output state vectors
field_names = fieldnames(workingstate) ;
for numfields = 1:length(field_names)
    eval([' state.' char( field_names(numfields) ) ' = workingstate.' char( field_names(numfields) ) '(trecords) ; '])
end

%%%%%% done message
% fprintf('Done: \n ')
% endtime = toc ;
% fprintf('time (s): %d \n', endtime )

% end

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

end

empty_rows = all(cellfun(@(x) isempty(x) || (isstring(x) && strlength(x)==0), results), 2);

results = results(~empty_rows, :);



%%  data processing

num_runs = size(results, 1);
time_lengths = zeros(num_runs, 1);

for i = 1:num_runs
    time_lengths(i) = length(results{i, 9}.time_myr); 
end

[~, max_idx] = max(time_lengths); 
base_time = results{max_idx, 9}.time_myr; 


interpolated_results = cell(num_runs, 1);


for i = 1:num_runs
    
    time_myr = results{i, 9}.time_myr; 
    state_data = results{i, 9}; 
    
    interpolated_state = struct();
    interpolated_state.time_myr = base_time; 


    fields = fieldnames(state_data);
    for f = 1:length(fields)
        field_name = fields{f};
        if ~isequal(field_name, 'time_myr') 
            field_data = state_data.(field_name);
            interpolated_state.(field_name) = interp1(time_myr, field_data, base_time, 'linear', NaN);
        end
    end

    interpolated_results{i} = interpolated_state;
end

d13c_DIC = NaN(num_runs, length(base_time));
T_p = NaN(num_runs, length(base_time));
Atmospheric_CO2_ppm = NaN(num_runs, length(base_time));
Atmospheric_CH4_a_ppm = NaN(num_runs, length(base_time));
f_volcanic_CO2 = NaN(num_runs, length(base_time));
f_CH4 = NaN(num_runs, length(base_time));
d13C_input = NaN(num_runs, length(base_time));

for i = 1:num_runs
    d13c_DIC(i, :) = 0.2 * (interpolated_results{i}.d13c_DIC_d + interpolated_results{i}.d13c_DIC_di + interpolated_results{i}.d13c_DIC_h + interpolated_results{i}.d13c_DIC_p + interpolated_results{i}.d13c_DIC_s);  
    T_p(i, :) = interpolated_results{i}.T_p;
    Atmospheric_CO2_ppm(i, :) = interpolated_results{i}.Atmospheric_CO2_ppm;
    Atmospheric_CH4_a_ppm(i, :) = interpolated_results{i}.Atmospheric_CH4_a_ppm;
    f_volcanic_CO2(i, :) = interpolated_results{i}.f_volcanic_CO2;
    f_CH4(i, :) = interpolated_results{i}.f_CH4;
    d13C_input(i, :) = interpolated_results{i}.d13c_CO2_input;
end

d13C_input_mean = nanmean(d13C_input, 1); 
d13C_input_5 = prctile(d13C_input, 5, 1); 
d13C_input_95 = prctile(d13C_input, 95, 1); 

f_volcanic_CO2_mean = nanmean(f_volcanic_CO2, 1); 
f_volcanic_CO2_5 = prctile(f_volcanic_CO2, 5, 1); 
f_volcanic_CO2_95 = prctile(f_volcanic_CO2, 95, 1); 

f_CH4_mean = nanmean(f_CH4, 1); 
f_CH4_5 = prctile(f_CH4, 5, 1); 
f_CH4_95 = prctile(f_CH4, 95, 1); 


d13c_DIC_mean = nanmean(d13c_DIC, 1); 
d13c_DIC_5 = prctile(d13c_DIC, 5, 1); 
d13c_DIC_95 = prctile(d13c_DIC, 95, 1); 

T_p_mean = nanmean(T_p, 1); 
T_p_5 = prctile(T_p, 5, 1); 
T_p_95 = prctile(T_p, 95, 1); 

Atmospheric_CO2_mean = nanmean(Atmospheric_CO2_ppm, 1); 
Atmospheric_CO2_5 = prctile(Atmospheric_CO2_ppm, 5, 1); 
Atmospheric_CO2_95 = prctile(Atmospheric_CO2_ppm, 95, 1); 

Atmospheric_CH4_mean = nanmean(Atmospheric_CH4_a_ppm, 1); 
Atmospheric_CH4_5 = prctile(Atmospheric_CH4_a_ppm, 5, 1); 
Atmospheric_CH4_95 = prctile(Atmospheric_CH4_a_ppm, 95, 1); 

plot_time = base_time * 1e6 ;

%%

color_line    = rand(1,3);
color_shade   = rand(1,3);
color_scatter = rand(1,3);

figure('Color', 'w', 'Position', [100 100 1200 1000]);
tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
ax1 = nexttile;
hold on;
box on;
a1 = plot(plot_time, d13c_DIC_mean,'-','LineWidth',2,'Color',color_line);
fill([(plot_time)', fliplr((plot_time)')],[(d13c_DIC_5), fliplr((d13c_DIC_95))],color_line,'LineStyle','none','facealpha',0.2);hold on;
fill([(d13C_p(:,1))', fliplr((d13C_p(:,1))')],[(d13C_p(:,2))', fliplr((d13C_p(:,4))')],color_scatter,'LineStyle','none','facealpha',0.2);hold on;
a2 = plot(d13C_p(:,1), d13C_p(:,3),'-','LineWidth',2,'Color',color_scatter);
legend([a1,a2],'CHEES','同化')
xlim([0 6000])
xlabel('Time (yr)', 'FontSize', 12)
ylabel(['δ^{13}C of DIC (', char(8240), ')'], 'FontSize', 12)
title('a', 'FontSize', 14, 'FontWeight', 'bold')

ax2 = nexttile;
hold on;
box on;
plot(plot_time, Atmospheric_CO2_mean,'-','LineWidth',2,'Color',color_line)
fill([(plot_time)', fliplr((plot_time)')],[(Atmospheric_CO2_5), fliplr((Atmospheric_CO2_95))],color_line,'LineStyle','none','facealpha',0.2);hold on;
fill([(pco2(:,1))', fliplr((pco2(:,1))')],[(pco2(:,2))', fliplr((pco2(:,4))')],color_scatter,'LineStyle','none','facealpha',0.2);hold on;
plot(pco2(:,1), pco2(:,3),'-','LineWidth',2,'Color',color_scatter)
xlim([0 6000])
xlabel('Time (yr)', 'FontSize', 12)
ylabel('Atm. CO_{2} (ppm)', 'FontSize', 12)
title('b', 'FontSize', 14, 'FontWeight', 'bold')

ax3 = nexttile;
hold on;
box on;
plot(plot_time, T_p_mean-273,'-','LineWidth',2,'Color',color_line)
fill([(plot_time)', fliplr((plot_time)')],[(T_p_5-273), fliplr((T_p_95-273))],color_line,'LineStyle','none','facealpha',0.2);hold on;
fill([(SST(:,1))', fliplr((SST(:,1))')],[(SST(:,2))', fliplr((SST(:,4))')],color_scatter,'LineStyle','none','facealpha',0.2);hold on;
plot(SST(:,1), SST(:,3),'-','LineWidth',2,'Color',color_scatter)
xlim([0 6000])
xlabel('Time (yr)', 'FontSize', 12)
ylabel('SST (°C)', 'FontSize', 12)
title('c', 'FontSize', 14, 'FontWeight', 'bold')

ax4 = nexttile;
hold on;
box on;
d1=plot(plot_time, f_CH4_mean,'-','LineWidth',2,'Color',color_line);
fill([(plot_time)', fliplr((plot_time)')],[(f_CH4_5), fliplr((f_CH4_95))],color_line,'LineStyle','none','facealpha',0.2);hold on;
d2=plot(plot_time, f_volcanic_CO2_mean,'-','LineWidth',2,'Color',color_shade);
fill([(plot_time)', fliplr((plot_time)')],[(f_volcanic_CO2_5), fliplr((f_volcanic_CO2_95))],color_shade,'LineStyle','none','facealpha',0.2);hold on;
xlim([0 6000])
legend([d1,d2],'CH_4','CO_2','box','off','location','northwest');
xlabel('Time (yr)', 'FontSize', 12)
ylabel('Carbon input  (mol/y)', 'FontSize', 12)
title('d', 'FontSize', 14, 'FontWeight', 'bold')

ax5 = nexttile;
hold on;
box on;
plot(plot_time, d13C_input_mean,'-','LineWidth',2,'Color',color_line);
fill([(plot_time)', fliplr((plot_time)')],[(d13C_input_5), fliplr((d13C_input_95))],color_line,'LineStyle','none','facealpha',0.2);hold on;
xlim([0 6000])
xlabel('Time (yr)', 'FontSize', 12)
ylabel(['δ^{13}C input (', char(8240), ')'], 'FontSize', 12)
title('e', 'FontSize', 14, 'FontWeight', 'bold')

ax6 = nexttile;
hold on;
box on;
plot(plot_time, Atmospheric_CH4_mean,'-','LineWidth',2,'Color',color_line);
fill([(plot_time)', fliplr((plot_time)')],[(Atmospheric_CH4_5), fliplr((Atmospheric_CH4_95))],color_line,'LineStyle','none','facealpha',0.2);hold on;
xlim([0 6000])
xlabel('Time (yr)', 'FontSize', 12)
ylabel('Atm. CH_{4} (ppm)', 'FontSize', 12)
title('f', 'FontSize', 14, 'FontWeight', 'bold')

%% === 导出用于画图的时间序列到 CSV ===
% 时间使用 plot_time（单位：year），和图上一致
% 只保留 t >= 0 的时间点
idx = plot_time >= 0;          % 逻辑索引，对所有序列统一使用

% 时间
time_yr = plot_time(idx);
time_yr = time_yr(:);          % 强制转成列向量

% SST（转为 ℃）
SST_mean_degC = (T_p_mean(idx) - 273);
SST_mean_degC = SST_mean_degC(:);

SST_5_degC = (T_p_5(idx) - 273);
SST_5_degC = SST_5_degC(:);

SST_95_degC = (T_p_95(idx) - 273);
SST_95_degC = SST_95_degC(:);

% δ13C_DIC
d13c_DIC_mean_col = d13c_DIC_mean(idx);
d13c_DIC_mean_col = d13c_DIC_mean_col(:);

d13c_DIC_5_col = d13c_DIC_5(idx);
d13c_DIC_5_col = d13c_DIC_5_col(:);

d13c_DIC_95_col = d13c_DIC_95(idx);
d13c_DIC_95_col = d13c_DIC_95_col(:);

% 大气 CO2
Atmospheric_CO2_mean_col = Atmospheric_CO2_mean(idx);
Atmospheric_CO2_mean_col = Atmospheric_CO2_mean_col(:);

Atmospheric_CO2_5_col = Atmospheric_CO2_5(idx);
Atmospheric_CO2_5_col = Atmospheric_CO2_5_col(:);

Atmospheric_CO2_95_col = Atmospheric_CO2_95(idx);
Atmospheric_CO2_95_col = Atmospheric_CO2_95_col(:);

% CH4 通量
f_CH4_mean_col = f_CH4_mean(idx);
f_CH4_mean_col = f_CH4_mean_col(:);

f_CH4_5_col = f_CH4_5(idx);
f_CH4_5_col = f_CH4_5_col(:);

f_CH4_95_col = f_CH4_95(idx);
f_CH4_95_col = f_CH4_95_col(:);

% 火山 CO2 通量
f_volcanic_CO2_mean_col = f_volcanic_CO2_mean(idx);
f_volcanic_CO2_mean_col = f_volcanic_CO2_mean_col(:);

f_volcanic_CO2_5_col = f_volcanic_CO2_5(idx);
f_volcanic_CO2_5_col = f_volcanic_CO2_5_col(:);

f_volcanic_CO2_95_col = f_volcanic_CO2_95(idx);
f_volcanic_CO2_95_col = f_volcanic_CO2_95_col(:);

% 输入 δ13C
d13C_input_mean_col = d13C_input_mean(idx);
d13C_input_mean_col = d13C_input_mean_col(:);

d13C_input_5_col = d13C_input_5(idx);
d13C_input_5_col = d13C_input_5_col(:);

d13C_input_95_col = d13C_input_95(idx);
d13C_input_95_col = d13C_input_95_col(:);

% 大气 CH4
Atmospheric_CH4_mean_col = Atmospheric_CH4_mean(idx);
Atmospheric_CH4_mean_col = Atmospheric_CH4_mean_col(:);

Atmospheric_CH4_5_col = Atmospheric_CH4_5(idx);
Atmospheric_CH4_5_col = Atmospheric_CH4_5_col(:);

Atmospheric_CH4_95_col = Atmospheric_CH4_95(idx);
Atmospheric_CH4_95_col = Atmospheric_CH4_95_col(:);

% 组装成一个 table（所有变量现在都是 N×1 列向量）
TS_out = table( ...
    time_yr, ...
    d13c_DIC_mean_col, d13c_DIC_5_col, d13c_DIC_95_col, ...
    Atmospheric_CO2_mean_col, Atmospheric_CO2_5_col, Atmospheric_CO2_95_col, ...
    SST_mean_degC, SST_5_degC, SST_95_degC, ...
    f_CH4_mean_col, f_CH4_5_col, f_CH4_95_col, ...
    f_volcanic_CO2_mean_col, f_volcanic_CO2_5_col, f_volcanic_CO2_95_col, ...
    d13C_input_mean_col, d13C_input_5_col, d13C_input_95_col, ...
    Atmospheric_CH4_mean_col, Atmospheric_CH4_5_col, Atmospheric_CH4_95_col, ...
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
% 写出到当前目录的 CSV 文件
writetable(TS_out, fullpath);
fprintf('  Saved %s\n', fullpath);
end
%% ---- Local functions (must be at end of script) ----
function x = truncnorm(mu, sigma, lo, hi)
x = mu + sigma*randn();
while x < lo || x > hi
    x = mu + sigma*randn();
end
end


