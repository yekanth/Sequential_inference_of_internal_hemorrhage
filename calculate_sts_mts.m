function [states,measurements,measurements_TN_FP] = calculate_sts_mts(VSUBJECTS_VALID,noise,noise_type,per_BV_critical,model_known,H,I,tran_time,filter_measurements,st,mt)
    
    n_VP = length(VSUBJECTS_VALID(:,1));
    
    if model_known == 1
        VSUBJECTS_VALID = hr_scale_parameters(VSUBJECTS_VALID);
        req_params = [VSUBJECTS_VALID(:,1) VSUBJECTS_VALID(:,2) VSUBJECTS_VALID(:,3) VSUBJECTS_VALID(:,11) VSUBJECTS_VALID(:,12)];
        BV = req_params(4);
        [states,measurements,measurements_TN_FP] = kinetics_transition(req_params,H,I,tran_time,per_BV_critical,noise,noise_type,filter_measurements,st); 
    end

    %Resampling Measurements
    if ~isempty(measurements)
        measurements.HCT = resample(measurements.HCT,0:st:length(measurements.HCT)*st-st,1/mt);
        measurements.y_true = resample(measurements.y_true,0:st:length(measurements.y_true)*st-st,1/mt);
        measurements.z_dot = resample(measurements.z_dot,0:st:length(measurements.z_dot)*st-st,1/mt);
    end
    
    if ~isempty(measurements_TN_FP)
        measurements_TN_FP.HCT = resample(measurements_TN_FP.HCT,0:st:length(measurements_TN_FP.HCT)*st-st,1/mt);
        measurements_TN_FP.y_true = resample(measurements_TN_FP.y_true,0:st:length(measurements_TN_FP.y_true)*st-st,1/mt);
        measurements_TN_FP.z_dot = resample(measurements_TN_FP.z_dot,0:st:length(measurements_TN_FP.z_dot)*st-st,1/mt);
    end


end