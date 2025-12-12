
function [new_d13C] = high_resolution_d13C_data(interp_time,d13C_origin)


% interp_time = [0.403: 0.0005 : 0.621];
d13C_age = d13C_origin(:,1);
d13C_data = d13C_origin(:,2);
% logical_indices_new_d13C_NCIE = new_d13C_age >=  0.491737 & new_d13C_age <= 0.516954;
logical_indices = d13C_age >=  min(interp_time) & d13C_age <= max(interp_time);
d13C_age_CIE = d13C_age(logical_indices);
d13C_data_CIE = d13C_data(logical_indices);
new_d13C_data = interp1(d13C_age_CIE, d13C_data_CIE, interp_time,'linear');
new_d13C_CIE = [interp_time',new_d13C_data'];
new_d13C = [d13C_origin;new_d13C_CIE];
new_d13C = sortrows(new_d13C);
[a_1, idx, ~] = unique(new_d13C(:, 1));
a_2 = new_d13C(idx, 2);
new_d13C=[a_1,a_2];
nan_indices = any(isnan(new_d13C), 2);  
new_d13C = new_d13C(~nan_indices, :);   


end