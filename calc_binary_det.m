function [IC,H_RMSE] = calc_binary_det(SN,H,I,trantime,bar_max,states_estimated,theta3_max,VSUBJECTS_NAT,measurements,nominal,algo_ind,noise,trail_window,BP_detection_threshold,BP_noise,noise_type,st,mt)
    
    % For each algorithm generate a vector in time with 1s and 0s, where 1
    % represent hemorrhage detected for the given time and 0 represent
    % hemorhage not detected for a given time
    
    H_est_avg = [];
    H_RMSE = [];
    
    
    if length(states_estimated.x_hat) == 0
        IC = zeros(floor(length(measurements.HCT)*mt/st),1);
    else 
        IC = zeros(length(states_estimated.x_hat(:,1)),1);
    end
    
    if any(any(algo_ind(:) == [1,2,3,4]))==1
        if length(states_estimated.x_hat(:,1))>length(bar_max.x1_bar)
            IC = zeros(length(bar_max.x1_bar(1,:)),1);
        else
            IC = zeros(length(states_estimated.x_hat(:,1)),1);
        end
    end
    
    if algo_ind == 1
       if ~isempty(find(states_estimated.x_hat(1:length(IC),1) > bar_max.x1_bar(1,1:length(IC))'))
            IC(find(states_estimated.x_hat(1:length(IC),1) > bar_max.x1_bar(1,1:length(IC))')) = 1;
       end
    end
    
    if algo_ind == 2
       if ~isempty(find(states_estimated.x_hat(1:length(IC),2) > bar_max.x2_bar(1,1:length(IC))'))
            IC(find(states_estimated.x_hat(1:length(IC),2) > bar_max.x2_bar(1,1:length(IC))')) = 1;
       end
    end
        
    if algo_ind == 3
       temp_state = states_estimated.x_hat(1:length(IC),3);
       if ~isempty(find(temp_state < -0.01))
            IC(find(temp_state < -0.01)) = 1;
            H_est_avg = movmean(-temp_state/nominal(:,2),[trail_window/st,0]);
            H_RMSE = sqrt(mean((H_est_avg(trantime/st:end)-H*ones(length(H_est_avg)-(trantime/st)+1,1)).^2))/H;
       end
    end

    if algo_ind == 4
       if ~isempty(find(states_estimated.x_hat(1:length(IC),4) > bar_max.x2_dot_bar(1,1:length(IC))'))
            IC(find(states_estimated.x_hat(1:length(IC),4) > bar_max.x2_dot_bar(1,1:length(IC))')) = 1;
       end
    end

    
    if algo_ind == 5
        H_est = run_EKF(SN,H,I,trantime,measurements,VSUBJECTS_NAT,nominal,noise,st,mt);
        H_est_avg = movmean(H_est,[trail_window/st,0]);
        fprintf('EKF State Estimation for Subject %d\n',SN)
        if ~isempty(find(H_est_avg > 0.005))
            IC(find(H_est_avg>0.005)) = 1;
        end
        if ~isempty(find(H_est_avg < 0))
            IC(find(H_est_avg < 0)) = 0;
        end
        H_RMSE = sqrt(mean((H_est_avg(trantime/st:end)-H*ones(length(H_est_avg)-(trantime/st)+1,1)).^2))/H;
    end
    
    if algo_ind == 6
        BP_est = estimate_BP(SN,H,I,trantime,measurements,VSUBJECTS_NAT,nominal,noise,st,mt);
        BP_est = add_noise(BP_est,BP_noise,noise_type);
        BP_est_avg = movmean(BP_est,[trail_window/st,0]);
        fprintf('Estimating BP for Subject %d\n',SN)
        if ~isempty(find(BP_est_avg < BP_detection_threshold*BP_est_avg(trantime/st)))
            IC(find(BP_est_avg<BP_detection_threshold*BP_est_avg(trantime/st))) = 1;
        end
    end
    
    if algo_ind == 7
        BP_est = estimate_BP(SN,H,I,trantime,measurements,VSUBJECTS_NAT,nominal,noise,st,mt);
        BP_est = add_noise(BP_est,BP_noise,noise_type);
        BP_est_avg = movmean(BP_est,[trail_window/st,0]);
        fprintf('Estimating BP for Subject %d\n',SN)
        if ~isempty(find(BP_est_avg < 0.8*BP_est_avg(trantime/st)))
            IC(find(BP_est_avg<0.8*BP_est_avg(trantime/st))) = 1;
        end
    end
    
    if algo_ind == 8
        BP_est = estimate_BP(SN,H,I,trantime,measurements,VSUBJECTS_NAT,nominal,noise,st,mt);
        BP_est = add_noise(BP_est,BP_noise,noise_type);
        BP_est_avg = movmean(BP_est,[trail_window/st,0]);
        fprintf('Estimating BP for Subject %d\n',SN)
        if ~isempty(find(BP_est_avg < 0.7*BP_est_avg(trantime/st)))
            IC(find(BP_est_avg<0.7*BP_est_avg(trantime/st))) = 1;
        end
    end
    
    if algo_ind == 9
        HCT_est = measurements.HCT;
        HCT_est_avg = movmean(HCT_est,[trail_window/mt,0]);
        fprintf('Running HCT Naive for Subject %d\n',SN)
        IC = zeros(ceil(length(measurements.HCT)*mt),1);
        if ~isempty(find(HCT_est_avg < 0.9*HCT_est_avg(floor(trantime/mt))))
            IC(floor(find(HCT_est_avg<0.9*HCT_est_avg(floor(trantime/mt))))) = 1;
        end
        IC = resample(IC,0:mt:length(IC)*mt-mt,1/st);
    end
    
    if algo_ind == 10
        HCT_est = measurements.HCT;
        HCT_est_avg = movmean(HCT_est,[trail_window/mt,0]);
        fprintf('Running HCT Naive for Subject %d\n',SN)
        IC = zeros(ceil(length(measurements.HCT)*mt),1);
        if ~isempty(find(HCT_est_avg < 0.8*HCT_est_avg(floor(trantime/mt))))
            IC(floor(find(HCT_est_avg<0.8*HCT_est_avg(floor(trantime/mt))))) = 1;
        end
        IC = resample(IC,0:mt:length(IC)*mt-mt,1/st);
    end
    
    if algo_ind == 11
        HCT_est = measurements.HCT;
        HCT_est_avg = movmean(HCT_est,[trail_window/mt,0]);
        fprintf('Running HCT Naive for Subject %d\n',SN)
        IC = zeros(ceil(length(measurements.HCT)*mt),1);
        if ~isempty(find(HCT_est_avg < 0.7*HCT_est_avg(floor(trantime/mt))))
            IC(floor(find(HCT_est_avg<0.7*HCT_est_avg(floor(trantime/mt))))) = 1;
        end
        IC = resample(IC,0:mt:length(IC)*mt-mt,1/st);
    end
    
    if algo_ind == 12
        H_est = run_UKF(SN,H,I,trantime,measurements,VSUBJECTS_NAT,nominal,noise,st,mt);
        H_est_avg = movmean(H_est,[trail_window/st,0]);
        fprintf('UKF State Estimation for Subject %d\n',SN)
        if ~isempty(find(H_est_avg > 0.005))
            IC(find(H_est_avg>0.005)) = 1;
        elseif ~isempty(find(H_est_avg < 0))
            IC(find(H_est_avg < 0)) = 0;
        end
        H_RMSE = sqrt(mean((H_est_avg(trantime/st:end)-H*ones(length(H_est_avg)-(trantime/st)+1,1)).^2))/H;
    end
    
    if algo_ind == 13
        %IC = zeros(length(measurements.HCT)*mt/st,1);
        theta3 = parameter_estimator(SN,H,I,trantime,theta3_max,measurements,VSUBJECTS_NAT,nominal,noise,st,mt);
        if ~isempty(find(theta3 > theta3_max))
            IC(find(theta3 > theta3_max)) = 1;
        end  
    end
        
end