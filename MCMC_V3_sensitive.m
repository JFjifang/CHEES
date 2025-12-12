
clc;clear


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%  load data and Data preprocessing %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%-------------------------load data-------------------------------------

global stepnumber
global pars
global workingstate

%%%% load d13C_p
d13C_p = load('data/d13Cp_data.txt');

%%%% load temperature
SST = load('data/SST.txt');

%%%% load pco2
pco2 = load('data/CO2.txt');

age = d13C_p(:,1);

 
step = 1 ;  
num_intervals = floor((length(age) - 1) / step);
segments_3 = zeros(num_intervals, 2);

for i = 1:num_intervals
    idx = (i - 1) * step + 1;  
    segments_3(i, 1) = age(idx);      
    segments_3(i, 2) = age(idx + step); 
end


if mod(length(age) - 1, step) ~= 0
    segments_3(end + 1, :) = [age(end - step), age(end)];
end

n_repeats = 500; 
segments = repmat(segments_3, n_repeats, 1);

num_iterations = size(segments,1); % Adjust the num  ber of iterations as needed
acceptance_rate = zeros(1, num_iterations);
proposed_accepted_state = {};
current_likelihoods = zeros( 1,num_iterations); 
proposed_likelihoods = zeros( 1,num_iterations); 
accepted_likelihoods = zeros( 1,num_iterations); 
max_consecutive_failures = 8;  % Set the maximum consecutive failures
consecutive_failures = 0;      % Initialize consecutive failures counter
perturbation_scaling_factor = 1;

emission_model = load('data/emission.txt');
[new_emission_model] = add_data_emission(emission_model,d13C_p);
time = new_emission_model(:,1)';
input_C = new_emission_model(:,2)';
d13C_input = new_emission_model(:,3)';
d13C_methane = new_emission_model(:,4)';

rng('shuffle');

d13C_CO2 = -5 ;     % -5 - 20 * rand(1);
% d13C_methane = -60 ; % -50 - 20 * rand(1);
methane_oxygen = 9e3 ; %+ 2.5e3 * rand(1);

k_CH4_diff = 0.7 + 0.2 * rand(1); % 0.1 - 0.5 生成后扩散到水柱中的甲烷
k_CH4_from_ocean = 0.8 + 0.2 * rand(1); % 0.2 - 0.8 有多少甲烷是来自海洋的
k_CH4_oxiwatercolum = 0.2 + 0.1 * rand(1); % 0.1 - 0.5 这些甲烷气泡在水柱中有多少被氧化


%% 

for iter = 1 : size(segments, 1)

        segment_start = segments(iter, 1);
        segment_end = segments(iter, 2);
        disp(['MCMC iteration: ', num2str(iter)]);
        [current_pars,state] = CHEES_frontend(time,input_C,d13C_input,d13C_CO2,d13C_methane,methane_oxygen,k_CH4_diff,k_CH4_from_ocean,k_CH4_oxiwatercolum); 
        [current_likelihood,ll_d13,ll_pco2,ll_sst] = log_likelihood(state.time_myr, d13C_p, state.d13c_DIC_p, pco2, state.Atmospheric_CO2_ppm, SST,state.GAST);              
        current_likelihoods(iter) = current_likelihood;
        new_input_C_array = []; 
        new_d13C_input_array = [];
        proposed_likelihood_array = [];
        
    parfor i = 1:28
%         disp(['proposed iteration:', num2str(i)]);
        [new_input_C, new_d13C_input,new_d13C_CH4] = propose_new_parameters(time,input_C,d13C_input,d13C_methane,perturbation_scaling_factor);
       % figure;plot(new_input_C)
        [~, proposed_state] = CHEES_frontend(time, new_input_C, new_d13C_input,d13C_CO2,new_d13C_CH4,methane_oxygen,k_CH4_diff,k_CH4_from_ocean,k_CH4_oxiwatercolum); 
        [proposed_likelihood,~,~,~] = log_likelihood(proposed_state.time_myr,d13C_p,proposed_state.d13c_DIC_p,pco2,proposed_state.Atmospheric_CO2_ppm,SST,proposed_state.GAST);
        new_input_C_array(i, :) = new_input_C;
        new_d13C_input_array(i, :) = new_d13C_input;
        new_d13C_CH4_array(i, :) = new_d13C_CH4;
        proposed_likelihood_array(i) = proposed_likelihood;
        proposed_state_cell{i}= proposed_state;
    end 

    
    % Find the maximum likelihood and corresponding parameters among the 500
    [max_proposed_likelihood, max_index] = max(proposed_likelihood_array);
    best_new_input_C = new_input_C_array(max_index, :);
    best_new_d13C_input = new_d13C_input_array(max_index, :);
    best_new_CH4_input = new_d13C_CH4_array(max_index, :);

    proposed_likelihoods(iter) = max_proposed_likelihood;

  if max_proposed_likelihood > current_likelihood || log(rand()) < (max_proposed_likelihood - current_likelihood)
        last_likelihood = max_proposed_likelihood;
        acceptance_rate(iter) = 1;  % Record acceptance rate
        input_C = best_new_input_C;
        d13C_input = best_new_d13C_input;
        d13C_methane = best_new_CH4_input;
        accepted_likelihoods = [accepted_likelihoods, max_proposed_likelihood];
        accepted_proposed_state = proposed_state_cell{max_index};
        proposed_accepted_state = [proposed_accepted_state; {accepted_proposed_state}];
        consecutive_failures = 0;  % Reset consecutive failures counter
    else
        % Reject the proposal
        acceptance_rate(iter) = 0;
        consecutive_failures = consecutive_failures + 1;
        % Check for consecutive rejections
        if consecutive_failures >= max_consecutive_failures
            perturbation_scaling_factor = 0.8 * perturbation_scaling_factor; % Reduce by 20%
            consecutive_failures = 0;  

        end      
  end
    
       
end


%% 

color_scatter = [0.284215890004489	0.784872858146586	0.380570272702924 ] ; 


total_carbon_emission = trapz(state.time_myr*1e6,state.CO2_input*12*10^-15)
total_CH4_emission    = trapz(state.time_myr*1e6,state.f_CH4*12*10^-15)
total_CO2_emission    = trapz(state.time_myr*1e6,state.f_volcanic_CO2*12*10^-15)


figure;
subplot(2,3,1)
hold on
box on
% stairs(tate.time_myr,state.CO2_input,'k');
a1=plot(state.time_myr,state.f_CH4,'-o');
a2=plot(state.time_myr,state.f_volcanic_CO2,'b');
legend([a1,a2],'methane','volcanic CO_2');
xlim([0 0.0063])
% ylim([0 62e12])
xlabel('Time (Ma)')
ylabel('CO_{2} input (mol/yr)')
title('a', 'FontSize', 14, 'FontWeight', 'bold')

subplot(2,3,2)
hold on
box on
plot(state.time_myr,state.d13c_CO2_input,'k');
xlim([0 0.0063])
xlabel('Time (Ma)', 'FontSize', 12)
ylabel(['δ^{13}C of CO_{2} input (', char(8240), ')'], 'FontSize', 12)
title('b', 'FontSize', 14, 'FontWeight', 'bold')


subplot(2,3,3)
hold on
box on
plot(state.time_myr,state.d13C_CH4,'k');
xlim([0 0.0063])
xlabel('Time (Ma)', 'FontSize', 12)
ylabel(['δ^{13}C of CO_{2} input (', char(8240), ')'], 'FontSize', 12)
title('c', 'FontSize', 14, 'FontWeight', 'bold')


%%%% CO2 (ppm)
subplot(2,3,4)
hold on
box on
plot(state.time_myr,state.Atmospheric_CO2_ppm,'r')
fill([(pco2(:,1)/1e6)', fliplr((pco2(:,1)/1e6)')],[(pco2(:,2))', fliplr((pco2(:,4))')],color_scatter,'LineStyle','none','facealpha',0.2);hold on;
plot(pco2(:,1)/1e6, pco2(:,3),'-','LineWidth',2,'Color',color_scatter)
xlim([0 0.0063])
xlabel('Time (Ma)', 'FontSize', 12)
ylabel('Atm. CO_{2} (ppm)', 'FontSize', 12)
title('d', 'FontSize', 14, 'FontWeight', 'bold')


subplot(235)
hold on
box on 
plot(state.time_myr, state.d13c_DIC_p,'r-')
fill([(d13C_p(:,1)/1e6)', fliplr((d13C_p(:,1)/1e6)')],[(d13C_p(:,2))', fliplr((d13C_p(:,4))')],color_scatter,'LineStyle','none','facealpha',0.2);hold on;
plot(d13C_p(:,1)/1e6, d13C_p(:,3),'-','LineWidth',2,'Color',color_scatter);
xlim([0 0.0063])
xlabel('Time (Ma)', 'FontSize', 12)
ylabel(['δ^{13}C of POC (', char(8240), ')'], 'FontSize', 12)
title('e', 'FontSize', 14, 'FontWeight', 'bold')

subplot(2,3,6)
hold on
box on
fill([(SST(:,1)/1e6)', fliplr((SST(:,1)/1e6)')],[(SST(:,2))', fliplr((SST(:,4))')],color_scatter,'LineStyle','none','facealpha',0.2);hold on;
plot(SST(:,1)/1e6, SST(:,3),'-','LineWidth',2,'Color',color_scatter)
plot(state.time_myr,state.T_p-273,'r-')
xlim([0 0.0063])
xlabel('Time (Ma)', 'FontSize', 12)
ylabel('T (^{0}C)', 'FontSize', 12)
title('f', 'FontSize', 14, 'FontWeight', 'bold')

