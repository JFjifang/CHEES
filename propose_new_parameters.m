% function [new_input_C, new_d13C_input,new_ccdeg] = propose_new_parameters(time,input_C,d13C_input,ccdeg,perturbation_scaling_factor,segment_start, segment_end,iter)
function [new_input_C, new_d13C_input,new_d13C_CH4] = propose_new_parameters(time,input_C,d13C_input,d13C_CH4,perturbation_scaling_factor)

        time_indices = time > 0  & time <= 6000 ;

        perturbation_input_C = 300e12; % 5e11
        perturbation_d13C_input = 40; % 10
        perturbation_d13C_CH4 = -30;

        proposed_params_input_C = input_C;
        proposed_params_d13C_input = d13C_input;
        proposed_params_d13C_CH4 = d13C_CH4;

        perturbation_randn = randn(1, 3);
        proposed_params_input_C(time_indices) = input_C(time_indices) + perturbation_input_C * perturbation_randn(1) * perturbation_scaling_factor;
        proposed_params_d13C_input(time_indices) = d13C_input(time_indices) + perturbation_d13C_input * perturbation_randn(2) * perturbation_scaling_factor;
        proposed_params_d13C_CH4(time_indices) = d13C_CH4(time_indices) + perturbation_d13C_CH4 * perturbation_randn(3)* perturbation_scaling_factor;

        d13C_CO2 = -5; 
        % d13C_methane = -60;
   
        for i = 1:length(proposed_params_d13C_input)
            proposed_params_d13C_input(i) = min(proposed_params_d13C_input(i), d13C_CO2);
        end
    
        % 对 proposed_params_d13C_input  限制
        proposed_params_d13C_CH4 = max(min(proposed_params_d13C_CH4, -30), -60);
        
         for i = 1:length(proposed_params_d13C_CH4)
            proposed_params_d13C_input(i) = max(proposed_params_d13C_input(i), proposed_params_d13C_CH4(i));
        end   

        proposed_params_input_C = max(proposed_params_input_C, 0);
    
        new_input_C = proposed_params_input_C;
        new_d13C_input = proposed_params_d13C_input;
        new_d13C_CH4 = proposed_params_d13C_CH4;

end
