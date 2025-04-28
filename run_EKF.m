function H_est = run_EKF(SN,H,I,trantime,measurements,VSUBJECTS_NAT,nominal,parameters_known,noise,st,mt)
    
    store_only_hemorrhage_state = 1;
    
    HCT_measured = measurements.HCT;
    theta = VSUBJECTS_NAT;
    
    alltime = length(HCT_measured)*mt;
    
    alpha_u_std = std(theta(:,1));
    alpha_h_std = std(theta(:,2));
    k_std = std(theta(:,3));
    V0_std = std(theta(:,11));
    H0_std = std(theta(:,12));
    
    theta_s = nominal;

    if parameters_known
        alpha_u = VSUBJECTS_NAT(SN,1);
        alpha_h = VSUBJECTS_NAT(SN,2);
        k = VSUBJECTS_NAT(SN,3);
        V0 = VSUBJECTS_NAT(SN,11);
        H0 = VSUBJECTS_NAT(SN,12);

        theta_s(:,1) = VSUBJECTS_NAT(SN,1);
        theta_s(:,2) = VSUBJECTS_NAT(SN,2);
        theta_s(:,3) = VSUBJECTS_NAT(SN,3);
        theta_s(:,11) = VSUBJECTS_NAT(SN,11);
        theta_s(:,12) = VSUBJECTS_NAT(SN,12);
    else
        alpha_u = theta_s(:,1);
        alpha_h = theta_s(:,2);
        k = theta_s(:,3);
        V0 = theta_s(:,11);
        H0 = theta_s(:,12);
    end
    
    %Get True States
    
    Irate = I.*[zeros(trantime/st,1); ones(round(alltime/st)-trantime/st,1)];
    Hrate = H.*[zeros(trantime/st,1); ones(round(alltime/st)-trantime/st,1)];
    H = [];
    
    ic = [0 0];
    
    %Uncmment if you need pts
    %[~,HCT_true,x_T] = NLF_hr_true_model(theta(SN,:),Hrate,Irate,st,ic);
    
    % EKF Code Starts
    
    
    %Covariance
    
    Q_P = [k_std^2 0 0;0 alpha_u_std^2 0;0 0 alpha_h_std^2];
    
    %State Covariance
    
    %P = 0.1*[0.001 0 0 0;0 0.001 0 0;0 0 0.001 0;0 0 0 0.001]; %State Covriance Matrix
    P = [0.01 0 0 0; 0 0.01 0 0;0 0 0.01 0;0 0 0 0.01];
    %Noise Covariance
    
    if noise==0
        R = 1e-4;
    else
        R = 0.1*noise;
    end
    
    JI = Irate; 
    
    x = [0,0,0,HCT_measured(1)*V0];
    
    HCT_hat(1) = HCT_measured(1);
    x_hat = x;
    h_state = [];
    
    for i = 2:length(Irate)
    
        %Prediction
        x_P = NLF_hr_transition_model(x,theta_s,JI(i-1),st);


    
        J = [-x(1) 0 0 0; (-x(3)/(1+alpha_h))+(JI(i-1)/(1+alpha_u)) -k*JI(i-1)/(1+alpha_u)^2 k*x(3)/(1+alpha_h)^2 0; 0 0 0.1 0];
        
        Q = J'*Q_P*J;
        
        A = [-k 1 -1 0; 0 0 -k/(1+alpha_h) 0; 0 0 0 0; (x(3)*x(4))/(x(1)+V0)^2 0 -x(4)/(x(1)+V0) -x(3)/(x(1)+V0)];
    
        P_P = P + st*[A*P + P*A' + Q];
    
        %Only store the required hemorrhage estimate
        if store_only_hemorrhage_state == 1
            h_state = [h_state;x_P(:,3)];
        else
            x_hat = [x_hat;x_P];
        end
        
        HCT_hat(i) = NLF_hr_measurement_model(x_P,theta_s);
    
    
    
        %Correction

        if mod(i,mt/st)==0
    
            H = [-x_P(4)/(x_P(1)+V0)^2 0 0 1/(x_P(1)+V0)];
    
            K = P_P*H'*(inv([H*P_P*H' + R]));
    
            x = x_P + K'.*(HCT_measured(floor(i*st/mt)) - HCT_hat(i));
        
            P = [eye(4) - K*H]*P_P*[eye(4) - H'*K'] + R*K*K';
        end
        
    end
    
    if store_only_hemorrhage_state == 1
        H_est = h_state;
    else
        H_est = x_hat(:,3);
    end
    
    % figure(2)
    % hold on;
    % plot(h_state)
    
    
    t = 0:st:alltime-st;
  
% figure(3)
% figure('units','normalized','outerposition',[0 0 1 1])
% subplot(2,3,1)
% plot(t,Irate)
% hold on;
% plot(t,-Hrate)
% xlabel('Time (min)')
% ylabel('Fluid Profile (L/min)')
% legend('Infusion','Hemorrhage')
% 
% subplot(2,3,2)
% plot(t,x_T(:,1))
% hold on;
% plot(t,x_hat(:,1))
% xlabel('Time (min)')
% ylabel('x1 true vs x1 estimated')
% legend('x1 true','x1 est')
% 
% subplot(2,3,3)
% plot(t,x_T(:,2))
% hold on;
% plot(t,x_hat(:,2))
% xlabel('Time (min)')
% ylabel('x2 true vs x2 estimated')
% legend('x2 true','x2 est')
% subplot(2,3,4)
% plot(t,Hrate)
% hold on;
% plot(t,x_hat(:,3))
% xlabel('Time (min)')
% ylabel('x3 true vs x3 estimated')
% legend('x3 true','x3 est')
% subplot(2,3,5)
% plot(t,x_T(:,3))
% hold on;
% plot(t,x_hat(:,4))
% ylabel('x4 true vs x4 estimated')
% legend('x4 true','x4 est')
% 
% subplot(2,3,6)
% plot(t,HCT_true)
% hold on;
% plot(t,HCT_hat)
% xlabel('Time (min)')
% ylabel('HCT true vs HCT estimated')
% legend('HCT true','HCT est')


    


    
end


