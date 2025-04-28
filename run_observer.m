function states_estimated = run_observer(vp_actual_measurements,nominal,n_VP,H,I,tran_time,min_sim_time,st,mt)
    
    % Function to estimate states using linear observer (Luenberger
    % observer)
    
    alpha_u = nominal(1);
    alpha_h = nominal(2);
    K = nominal(3);
    BV0 = nominal(11);
    HCT0 = nominal(12);
    
    %Use Butterworth poles for the Observer
    fc = 300;
    fs = 1000;
    [~,ppp,~] = butter(3,fc/(fs/2));
    
    %Plant Dynamics+
    A = [-K 1 -1; 0 0 -K/(1+alpha_h); 0 0 0 ];
    C = [1/BV0 0 0];
    
    %Find Observer Gains
    L =  place(A',C',[ppp(1) ppp(2) ppp(3)]);
    l1 = L(1);
    l2 = L(2);
    l3 = L(3);
    
    %For each VP, find the observed states and store into a string
    for i = 1:n_VP
        field = sprintf('subject_%d', i);
        this_VP_time = length(vp_actual_measurements.(field).y_true)*mt;
        y = vp_actual_measurements.(field).y_true;
        %Sanity Check
%         if length(y) ~= floor(this_VP_time/mt)
%             fprintf('Error - Stop Simulation - Lengths not equal in Observer')
%         end
        fprintf('Running Observer for Subject %d\n',i)
        x1 = 0;
        x2 = 0;
        x3 = 0;
        x_hat = [x1 x2 x3 0];
        x1_hat = 0;
        x2_hat = 0;
        x3_hat = 0;
        x1_dot_hat = 0;
        for j = 2:this_VP_time/st
            if j<tran_time/st
                inf = 0;
            else
                inf = I;
            end

            if mod(j,mt/st)==0
                x1 = x1 + st* (((-K - l1/BV0)*x1 + x2 - x3) + inf + l1*y(floor(j*st/mt)));
                x2 = x2 + st* (-l2*x1/BV0 - (K/(1+alpha_h))*x3 + (K*inf/(1+alpha_u)) + l2*y(floor(j*st/mt))); 
                x3 = x3 + st* (-l3*x1/BV0 + l3*y(floor(j*st/mt)));
                x1_hat = [x1_hat;x1];
                x2_hat = [x2_hat;x2];
                x3_hat = [x3_hat;x3];
                x1_dot_hat = [x1_dot_hat;((-K - l1/BV0)*x1 + x2 - x3) + inf + l1*y(floor(j*st/mt))];
            else
                x1 = x1 + st* (((-K)*x1 + x2 - x3) + inf);
                x2 = x2 + st* ( - (K/(1+alpha_h))*x3 + (K*inf/(1+alpha_u))); 
                x3 = x3 + st* (0);
                x1_hat = [x1_hat;x1];
                x2_hat = [x2_hat;x2];
                x3_hat = [x3_hat;x3];
                x1_dot_hat = [x1_dot_hat;((-K)*x1 + x2 - x3) + inf];
            end
            j;
        end
        x_hat = [x1_hat x2_hat x3_hat x1_dot_hat];
        states_estimated.(field).x_hat = x_hat;
        
    end
    
    %Plot for Sanity Check
    % figure(2)
    % subplot(4,1,1)
    % for i=1:n_VP
    %     field = sprintf('subject_%d', i);
    %     plot(states_estimated.(field).x_hat(1:floor(min_sim_time/st),1))
    %     hold on;
    % end
    % 
    % subplot(4,1,2)
    % for i=1:n_VP
    %     field = sprintf('subject_%d', i);
    %     plot(states_estimated.(field).x_hat(1:floor(min_sim_time/st),2))
    %     hold on;
    % end
    % 
    % subplot(4,1,3)
    % for i=1:n_VP
    %     field = sprintf('subject_%d', i);
    %     plot(states_estimated.(field).x_hat(1:floor(min_sim_time/st),3))
    %     hold on;
    % end
    % 
    % subplot(4,1,4)
    % for i=1:n_VP
    %     field = sprintf('subject_%d', i);
    %     plot(states_estimated.(field).x_hat(1:floor(min_sim_time/st),4))
    %     hold on;
    % end
             

end