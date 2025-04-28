function [ H ] = NLF_hr_measurement_model(x, theta_s)

    % Non-Linear Measurement model HCT measurement
    x1 = x(1); % Change in BV
    x4 = x(4); % RBC volume
     
    
    alpha_u = theta_s(:,1);
    alpha_h = theta_s(:,2);
    k      = theta_s(:,3);
    V0      = theta_s(:,11);
    H0      = theta_s(:,12);
    
    H = x4/(x1+V0); % HCT
    

end