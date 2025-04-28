function [TN,FP] = find_TN_FP(I,VSUBJECTS_VALID,theta3_max,algo_ind,tran_time,min_sim_time,per_BV_critical,model_known,noise,noise_type,window,prob_threshold,bar_max,nominal,trail_window,BP_detection_threshold,BP_noise,st,mt)
        % Function to calculate TN and FP based on system states and state
        % estimations from various algorithms 

        % First, normalized detection time is calculated and then using it
        % to differentiate TN, and FP
        H = 0;
        n_VP = length(VSUBJECTS_VALID(:,1));
        
        for i = 1:n_VP
                [~,~,measurements] = calculate_sts_mts(VSUBJECTS_VALID(i,:),noise,noise_type,per_BV_critical,model_known,H,I,tran_time,st,mt);
                field = sprintf('subject_%d', i);
                vp_actual_measurements.(field) = measurements;
                fprintf('Compiling Measurements for TN_FP Subject %d\n',i)
        end
                
        if any(any(algo_ind(:) == [1,2,3,4])) == 1
            states_estimated = run_observer(vp_actual_measurements,nominal,n_VP,H,I,tran_time,min_sim_time,st,mt);
        else
            for i = 1:n_VP
                field = sprintf('subject_%d', i);
                states_estimated.(field).x_hat = [];
            end
        end
        
        for i = 1:n_VP
            field = sprintf('subject_%d', i);
            VSUBJECTS_NAT = hr_scale_parameters(VSUBJECTS_VALID);
            req_params = [VSUBJECTS_NAT(i,1) VSUBJECTS_NAT(i,2) VSUBJECTS_NAT(i,3) VSUBJECTS_NAT(i,11) VSUBJECTS_NAT(i,12)];
            [~, ~, temp_measures] = kinetics_transition(req_params,H,I,tran_time,per_BV_critical,noise,noise_type,st);
            if ~isempty(temp_measures)
                temp_measures.HCT = resample(temp_measures.HCT,0:st:length(temp_measures.HCT)*st-st,1/mt);
                temp_measures.y_true = resample(temp_measures.y_true,0:st:length(temp_measures.y_true)*st-st,1/mt);
                temp_measures.z_dot = resample(temp_measures.z_dot,0:st:length(temp_measures.z_dot)*st-st,1/mt);
            end
            T_critical(i) = length(temp_measures.y_true)*mt;
            fprintf('Calculating TN and FP for Subject %d\n',i)    
            for j = 1:length(algo_ind)
                field2 = sprintf('algo_%d', algo_ind(j));
                [IC.(field2),~] = calc_binary_det(i,H,I,tran_time,bar_max,states_estimated.(field),theta3_max,VSUBJECTS_NAT,temp_measures,nominal,algo_ind(j),noise,trail_window,BP_detection_threshold,BP_noise,noise_type,st,mt);
                %Using the binary detection array, calculate NDT
                [NDT.(field).(field2),DT.(field).(field2)] = calculate_NDT(IC.(field2),tran_time,window,prob_threshold,T_critical(i),st);
            end
        end
        
        for j = 1:length(algo_ind)
            FP_temp = 0;
            NDT_temp = [];
            field2 = sprintf('algo_%d', algo_ind(j));
            for i = 1:n_VP
                field = sprintf('subject_%d', i);
                FP_temp = FP_temp + length((NDT.(field).(field2)<1));
                if length(NDT.(field).(field2))==0
                    NDT_temp(i) = NaN;
                else
                    NDT_temp(i) = NDT.(field).(field2);
                end
                if length(DT.(field).(field2))==0
                    DT_temp(i) = NaN;
                else
                    DT_temp(i) = DT.(field).(field2);
                end
            end
            FP.(field2) = FP_temp;
            TN.(field2) = n_VP-FP_temp;
            NDT.(field2) = NDT_temp;
            DT.(field2) = DT_temp;
        end
    
end