function H_est = run_UKF(SN,H,I,trantime,measurements,VSUBJECTS_NAT,nominal,noise,st,mt)
    
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
    
    alpha_u = theta_s(:,1);
    alpha_h = theta_s(:,2);
    k = theta_s(:,3);
    V0 = theta_s(:,11);
    H0 = theta_s(:,12);
    
    %Get True States
    
    Irate = I.*[zeros(trantime/st,1); ones(round(alltime/st)-trantime/st,1)];
    Hrate = H.*[zeros(trantime/st,1); ones(round(alltime/st)-trantime/st,1)];
    H = [];
    
    ic = [0 0];
    
    %Uncomment if you need pts
    %[~,HCT_true,x_T] = NLF_hr_true_model(theta(SN,:),Hrate,Irate,st,ic);
    
    % EKF Code Starts
    
    
    %Covariance
    
    Q = 0.0000000001*eye(4);
    
    %State Covariance
    
    %P = 0.1*[0.001 0 0 0;0 0.001 0 0;0 0 0.001 0;0 0 0 0.001]; %State Covriance Matrix
    P = 10*[0.001 0 0 0;0 0.001 0 0;0 0 0.001 0;0 0 0 0.001];
    
    %Noise Covariance
    R = 0.1*noise;
    
    JI = Irate; 
    
    x = [0,0,0,HCT_measured(1)*V0]';
    
    HCT_hat(1) = HCT_measured(1);
    x_hat = x';
    h_state = [];
    
    % Define the sigma points and weights
    L=numel(x);                                 %numer of states
    alpha=1e-2;                                 %default, tunable
    ki=-1;                                       %default, tunable
    beta=5;                                     %default, tunable
    lambda=alpha^2*(L+ki)-L;                    %scaling factor
    c=L+lambda;                                 %scaling factor
    Wm=[lambda/c 0.5/c+zeros(1,2*L)];           %weights for means
    Wc=Wm;
    Wc(1)=Wc(1)+(1-alpha^2+beta);               %weights for covariance
    c=sqrt(c);

    % UKF State Estimation

    for i = 2:length(Irate)
        
        X = sigmas(x, P, c);
        
        x_bar = zeros(L, 1);
            
        for j = 1:2*L+1
            x_bar = x_bar + Wm(j)*NLF_hr_transition_model(X(:, j),theta_s,JI(i-1),st)';
        end
        
        P_bar = zeros(L, L);
        for j = 1:2*L+1
            P_bar = P_bar + Wc(j)*(NLF_hr_transition_model(X(:, j),theta_s,JI(i-1),st)' - x_bar)*(NLF_hr_transition_model(X(:, j),theta_s,JI(i-1),st)' - x_bar)';
        end
        
        P_bar = P_bar + Q;
        
        % Measurement update

        if mod(i,mt/st)==0
            z_bar = 0;
            for j = 1:2*L+1
                z_bar = z_bar + Wm(j)*NLF_hr_measurement_model(X(:, j),theta_s);
            end
        
            HCT_hat(i) = z_bar;
            
            P_zz = 0;
            for j = 1:2*L+1
                P_zz = P_zz + Wc(j)*(NLF_hr_measurement_model(X(:, j),theta_s) - z_bar)*(NLF_hr_measurement_model(X(:, j),theta_s) - z_bar)';
            end
            P_zz = P_zz + R;
        
            P_xz = zeros(L, 1);
            for j = 1:2*L+1
                P_xz = P_xz + Wc(j)*(NLF_hr_transition_model(X(:, j),theta_s,JI(i-1),st)' - x_bar)*(NLF_hr_measurement_model(X(:, j),theta_s) - z_bar)';
            end
        
            K = P_xz*inv(P_zz);
            x = x_bar + K*(HCT_measured(floor(i*st/mt)) - z_bar);
            P = P_bar - K*P_zz*K';
        end
    
        
        %Only store the required hemorrhage estimate
        if store_only_hemorrhage_state == 1
            h_state = [h_state;x(3)];
        else
            x_hat = [x_hat;x'];
        end
        
        
    end
    
    if store_only_hemorrhage_state == 1
        H_est = h_state;
    else
        H_est = x_hat(:,3);
    end
 

    
end

% Function to sample sigma points for UKF

function X = sigmas(x, P, c)
    L = length(x);
    X = zeros(L, 2*L+1);
    X(:, 1) = x;
    A = c*chol(P)';
    %A = sqrt(c*P);
    for i = 1:L
        X(:, i + 1) = x + A(:, i);
        X(:, i + L + 1) = x - A(:, i);
    end
end


