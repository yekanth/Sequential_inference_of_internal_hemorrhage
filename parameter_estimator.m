function theta3 = parameter_estimator(SN,H,I,trantime,theta3_max,measurements,VSUBJECTS_NAT,nominal,noise,st,mt)
    
    % Function to estimate parameter by quantifying parameter uncertainity
    HCT_measured = measurements.HCT;
    z_dot = measurements.z_dot;
    v = measurements.v;
    z = measurements.y_true;

    theta_s = nominal;
    
    alpha_u_n = theta_s(:,1);
    alpha_h_n = theta_s(:,2);
    k_n = theta_s(:,3);
    BV0_n = theta_s(:,11);
    H0_n = theta_s(:,12);

    sim_time = length(HCT_measured)*mt;
    
    %Get True States
    
    Irate = I.*[zeros(trantime/st,1); ones(round(sim_time/st)-trantime/st,1)];
    Hrate = H.*[zeros(trantime/st,1); ones(round(sim_time/st)-trantime/st,1)];

    %Parameter estimation
    phi = zeros(3,length(z));
    for k = 1:length(z)
        phi(1,k) = -z(k);
        phi(2,k) = Irate(k);
        phi(3,k) = v(k);
    end

    theta = zeros(3,length(z));
    theta(:,1) = [k_n 1/BV0_n (1/BV0_n)*(k_n/(1+alpha_u_n))]';
    gamma = [2 0 0;0 2 0;0 0 2]';

    for l = 2:sim_time/st

        if mod(l,mt/st)==0
            theta(:,l) = theta(:,l-1) + 0.0001*gamma.*[z_dot(floor(l*st/mt)) - theta(:,l-1)'*phi(:,floor(l*st/mt))]*phi(:,floor(l*st/mt));
        else
            theta(:,l) = theta(:,l-1);
        end

        if theta(3,l) > 10
           break
        end
    end

    theta3 = theta(3,:)';


end