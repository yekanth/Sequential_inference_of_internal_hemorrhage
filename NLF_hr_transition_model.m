function next_state = NLF_hr_transition_model(x,theta_s,input,st)
        
        % State transition model, given previous state of the system and
        % model parameters
        Dt = st;
        x1 = x(1);
        x2 = x(2);
        x3 = x(3);
        x4 = x(4);
    
        JI = input;
    
        alpha_u = theta_s(:,1);
        alpha_h = theta_s(:,2);
        k  = theta_s(:,3);
        V0 = theta_s(:,11);
        H0  = theta_s(:,12);
    
    
    
        x1_next = x1 + Dt*(-k*x1 + x2 - x3 + JI);
        x2_next = x2 + Dt*((-k/(1+alpha_h))*x3 + (k/(1+alpha_u)*JI));
        x3_next = x3;
        x4_next = x4+ Dt*(-x3*x4/(x1+V0));
  
        next_state = [x1_next,x2_next,x3_next,x4_next];
end