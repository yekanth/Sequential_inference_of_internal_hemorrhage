hold on;
close all;
clear all;
clc;

% Give Deatils here for Hemorrhage Testing

include_only_valid = 0; %Choose only to include physiologically valid patients
trail_window = 2; %in min, window for moving average of measurements

%HCT
noise = 0.01; %Measurement noise in fraction
noise_type = 0; % 0 for Guassian, 1 for uniform, 2 for exponential

%BP
BP_noise = 3; %in mmHg
BP_detection_threshold = 0.9;

I_H = 0.1:0.2:2; % Range of Infusion rate/Hemorrhage rate to be tested
per_BV_critical = 0.25; % Percentage of Critical Blood Volume

mt = 0.01; % Measurement Frequency 
st = 0.01; % Simulation Frequency



window = 10; % In minutes, to calculate probability

prob_threshold = 0.50; % Probability threshold to detect Hemorrhage

% Pick the algorithms to be tested

% 1,2,3,4 - Linear Algorithms, 5 - Observer Based, 6,7 - EKF, UKF
algo_ind = [1,2,3,4,5,6,7,8,9,10,11,12,13]; 
%algo_ind = [3,5,12];
% Choose specific algorithms to be included in ensemble
ensemble_algo_ind = algo_ind;

% Is the model assumed to be known (Choose 1)
model_known = 1;

% Load Virtual Patients

load('VP_Data/VSUBJECTS_COLLOID_2sigma.mat')

% Load Nominal Patient (or the identified posterior mean from experimental data)

nominal = load('VP_Data/NOMINAL_COLLOID.mat').PHI_MU;
nominal = hr_scale_parameters(nominal);

% Eliminate Physiologically Invalid VPs
min_HC = 0;
max_HC = 0.5;
min_BP = 0;
max_BP = 100;
min_CO = 0;
max_CO = 10;


% Simulate the VPs to store their measurements (HCT, Cardiac Output, Blood
% Pressure)
i = 1;

input_data.Infusion_inp = 0.01*ones(18000,1);
input_data.Hemorrhage_inp = zeros(18000,1);
input_data.Urine_inp = zeros(18000,1);

% If user choses to include only physiologically valid VPs
if include_only_valid == 1
    for vsubject_no = 1:200
    params = VSUBJECTS(vsubject_no,:);
    outputs = run_model(params, input_data);
    S_HC = outputs.HC_estimated; % HCT
    S_CO = outputs.CO_estimated; % Cardiac Output
    S_BP = outputs.BP_estimated; % Blood Pressure
    
    % Check if the simulated VPs are in physiological range
    if isempty(find(S_HC<min_HC))&&isempty(find(S_HC>max_HC))&&isempty(find(S_BP<min_BP))&&isempty(find(S_BP>max_BP))&&isempty(find(S_CO<min_CO))&&isempty(find(S_CO>max_CO))
        VSUBJECTS_VALID(i,:) = VSUBJECTS(vsubject_no,:);
        i = i+1;
    end
    end
else
    VSUBJECTS_VALID = VSUBJECTS;
end

% Chose a subset of valid VPs (100 in this case)
VSUBJECTS_VALID = VSUBJECTS_VALID(1:2:200,:);


temp_VSUBJECTS_VALID = hr_scale_parameters(VSUBJECTS_VALID);

%For Algo 13, Param Uncertainity
%Calculating maximum theta 3 actual (theta 3 is the hemorrhage state
%estimated using parameter uncertainity algorithm)
theta3_act = [];
for i = 1:1:length(temp_VSUBJECTS_VALID(:,1))
 
        alpha_u = temp_VSUBJECTS_VALID(i,1);
        BV0 = temp_VSUBJECTS_VALID(i,11);
        aplha_h = temp_VSUBJECTS_VALID(i,2);
        k = temp_VSUBJECTS_VALID(i,3);
        
        theta3_act = [theta3_act (1/BV0)*(k/(1+alpha_u))];
end

[aa,indices]=sort(theta3_act,'descend');
maxnIndexes = indices(:,ceil(0.1*length(theta3_act)));
theta3_max = max(theta3_act(maxnIndexes));

alpha_h_min = min(temp_VSUBJECTS_VALID(:,2));

n_VP = length(VSUBJECTS_VALID(:,1)); % Number of VPs

%Give Transit Times (initial time period when there is no hemorrhage)
tran_time = 100; % in min

%0.001:0.003:0.05-0.003
h_ind = 1; % Index for hemorrhage rate

for H = [0.01:0.01:0.1] % Hemorrhage rates in L/min
    i_ind = 1; % Index for infusion rate
    
    for I = [0.1:0.2:2].*H % Infusion rates in L/min
        
        if I/H>0 % Run for cases only when I/H > 0
            bar_max = [];
            
            if any(any(algo_ind(:) == [1,2,3,4])) == 1
                %Simulate States without Hemorrhage
                %Simulate Actual Measurements & States for each VP
                for i = 1:n_VP
              
                    [states,measurements,~] = calculate_sts_mts(VSUBJECTS_VALID(i,:),noise,noise_type,per_BV_critical,model_known,H,I,tran_time,st,mt);
                    [bar_states,~,~] = calculate_sts_mts(VSUBJECTS_VALID(i,:),noise,noise_type,per_BV_critical,model_known,0,I,tran_time,st,mt);
                    field = sprintf('subject_%d', i);
                    vp_actual_states.(field) = states;
                    if length(bar_states.x1_bar)*st<length(measurements.HCT)*mt
                        measurements.HCT = measurements.HCT(1:floor(length(bar_states.x1_bar)*st/mt),:);
                        measurements.y_true = measurements.y_true(1:floor(length(bar_states.x1_bar)*st/mt),:);
                    end
                    
                    vp_actual_measurements.(field) = measurements;
                    vp_bar_states.(field) = bar_states;
                    fprintf('Compiling States and Measurements for Subject %d\n',i)
                end
                
                [bar_max,min_sim_time,max_sim_time] = calculate_bar_max(vp_actual_measurements,vp_bar_states,n_VP,st,mt);
                %Run the Observer for scenarios when there is no hemorrhage
            
                states_estimated = run_observer(vp_actual_measurements,nominal,n_VP,H,I,tran_time,min_sim_time,st,mt);
            
            end
            
            % Based on the user selection, run the suite of algorithms
            % (same code as above but measurement equation is modified for
            % algo_ind > 4
            

            if (any(any(algo_ind(:) == [5,6,7,8,9,10,11,12,13]))==1)&&(~any(any(algo_ind(:) == [1,2,3,4])) == 1)
                %Simulate States without Hemorrhage
                %Simulate Actual Measurements & States for each VP
                for i = 1:n_VP
                    [states,measurements,~] = calculate_sts_mts(VSUBJECTS_VALID(i,:),noise,noise_type,per_BV_critical,model_known,H,I,tran_time,st,mt);
                    field = sprintf('subject_%d', i);
                    vp_actual_states.(field) = states;
                    vp_actual_measurements.(field) = measurements;
                    states_estimated.(field).x_hat = [];
                    fprintf('Compiling States and Measurements for Subject %d\n',i)
                end
                
                max_sim_time = 0;
                min_sim_time = 1e8;
                for i=1:n_VP
                    field = sprintf('subject_%d', i);
                    sim_time = length(vp_actual_measurements.(field).HCT)*mt;
                    if sim_time>max_sim_time
                        max_sim_time = sim_time;
                    end
                    if sim_time<min_sim_time
                        min_sim_time = sim_time;
                    end
                 end
                
                
            end

            
            
            %Find TN and FP ie run without hemorrhage
            
            [TN,FP] = find_TN_FP(I,VSUBJECTS_VALID,theta3_max,algo_ind,tran_time,min_sim_time,per_BV_critical,model_known,noise,noise_type,window,prob_threshold,bar_max,nominal,trail_window,BP_detection_threshold,BP_noise,st,mt);
            
            %For each time step, simultaneously evaluate all algorithms and
            %their Ensemble
            
            %Calculate the binary detection for each condition and Critical
            %Time
            
            for i = 1:n_VP
                field = sprintf('subject_%d', i);
                VSUBJECTS_NAT = hr_scale_parameters(VSUBJECTS_VALID);
                req_params = [VSUBJECTS_NAT(i,1) VSUBJECTS_NAT(i,2) VSUBJECTS_NAT(i,3) VSUBJECTS_NAT(i,11) VSUBJECTS_NAT(i,12)];
                [~, temp_measures,~] = kinetics_transition(req_params,H,I,tran_time,per_BV_critical,noise,noise_type,st);

                if ~isempty(temp_measures)
                    temp_measures.HCT = resample(temp_measures.HCT,0:st:length(temp_measures.HCT)*st-st,1/mt);
                    temp_measures.y_true = resample(temp_measures.y_true,0:st:length(temp_measures.y_true)*st-st,1/mt);
                    temp_measures.z_dot = resample(temp_measures.z_dot,0:st:length(temp_measures.z_dot)*st-st,1/mt);
                end
                T_critical(i) = length(temp_measures.y_true)*mt;
                
                for j = 1:length(algo_ind)
                    field2 = sprintf('algo_%d', algo_ind(j));
                    [IC.(field2),H_NRMSE.(field2).(field)] = calc_binary_det(i,H,I,tran_time,bar_max,states_estimated.(field),theta3_max,VSUBJECTS_NAT,temp_measures,nominal,algo_ind(j),noise,trail_window,BP_detection_threshold,BP_noise,noise_type,st,mt);
                    %Using the binary detection array, calculate NDT
                    [NDT.(field).(field2),DT.(field).(field2)] = calculate_NDT(IC.(field2),tran_time,window,prob_threshold,T_critical(i),st);
                end
            end
            
            
            %Calculate TP and FN for each algorithm
       
            for j = 1:length(algo_ind)
                TP_temp = 0;
                NDT_temp = [];
                field2 = sprintf('algo_%d', algo_ind(j));
                for i = 1:n_VP
                    field = sprintf('subject_%d', i);
                    TP_temp = TP_temp + length(find(NDT.(field).(field2)<1));
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
                    if length(H_NRMSE.(field2).(field))==0
                        H_NRMSE_temp(i) = NaN;
                    else
                        H_NRMSE_temp(i) = H_NRMSE.(field2).(field);
                    end
                    
                end
                TP.(field2) = TP_temp;
                FN.(field2) = n_VP-TP_temp;
                NDT.(field2) = NDT_temp;
                DT.(field2) = DT_temp;
                H_NRMSE_alg.(field2) = H_NRMSE_temp;
                
            end
            
            % Compile the TP, FN, TN, FP count for each I-H pair and for all 100 virtual patients in a matrix 
            for j = 1:length(algo_ind)
                field = sprintf('algo_%d', algo_ind(j));
                TP_comp.(field).compilation(h_ind,i_ind) = TP.(field);
                FN_comp.(field).compilation(h_ind,i_ind) = FN.(field);
                TN_comp.(field).compilation(h_ind,i_ind) = TN.(field);
                FP_comp.(field).compilation(h_ind,i_ind) = FP.(field);
                NDT_comp.(field).compilation(h_ind,i_ind) = mean(NDT.(field),'omitnan');
                DT_comp.(field).compilation(h_ind,i_ind) = mean(DT.(field),'omitnan');
                if ~isempty(H_NRMSE_alg.(field))
                    H_NRMSE_comp.(field).compilation(h_ind,i_ind) = mean(H_NRMSE_alg.(field));
                end
            end
         
            
                
        end
        i_ind = i_ind+1;
    end
    
        h_ind = h_ind + 1
    
end
                            
                        
                        
                    
                
                