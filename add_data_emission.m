function [new_simulated_data] = add_data_emission(emission,d13C_p)

time = emission(:,1);
input_C = emission(:,2);
d13C = emission(:,3);
d13C_CH4 = emission(:,4);
d13C_CO2 = emission(:,5);

simulated_data = [time,input_C]; % 模拟数据
observed_data_d13C_P = d13C_p; % 实测数据
% observed_data_pco2 = pco2;
% observed_data_SST = SST;

% 获取模拟数据和实测数据的时间点
simulated_time = simulated_data(:, 1);
observed_time_d13C_P = observed_data_d13C_P(:, 1)*1e6;
% observed_time_pco2 = observed_data_pco2(:,1)*1e6;
% observed_time_SST =observed_data_SST(:,1)*1e6;

combined_time = [simulated_time; observed_time_d13C_P];
combined_time = combined_time(combined_time >= -7e6 & combined_time <= 0.33e6);
% combined_time = combined_time(combined_time >= -185e6  & combined_time <= -181e6);
combined_time = unique(combined_time); % 去重复
combined_time = sort(combined_time); % 排序
% combined_time=combined_time;

new_input_C = interp1(time, input_C, combined_time,'linear');
new_d13C = interp1(time, d13C, combined_time,'linear');
new_d13C_CH4 = interp1(time, d13C_CH4, combined_time,'linear');
new_d13C_CO2 = interp1(time, d13C_CO2, combined_time,'linear');

new_simulated_data = [combined_time,new_input_C,new_d13C,new_d13C_CH4,new_d13C_CO2];


end
