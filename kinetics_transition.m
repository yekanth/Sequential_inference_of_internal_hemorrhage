function [states,measurements,measurements_TN_FP] = kinetics_transition(req_params,H,I,trantime,per_BV_critical,noise,noise_type,st)
    
    
    filter_measurements = 1;
    
    alpha_u = req_params(1);
    alpha_h = req_params(2);
    K = req_params(3);
    BV0 = req_params(4);
    HCT0 = req_params(5);   
    
    x1 = 0;
    x2 = 0;
    
    u = I;
    h = 0;
    
    x1_bar = x1;
    x2_bar = x2;
    x1_dot_bar = -K*x1 + x2 + u - h;
    x2_dot_bar = K*u/(1+alpha_u) - K*h/(1+alpha_h);
    HCT_bar = HCT0;
    
    %For Parameter method
    v = 0; %Integration of u
    v_bar = 0;
    z_dot = 0; %Differentiation of fraction of BV from HCT
    z = 0;
    
    Vr = HCT0*BV0;
    Vr_temp = 0;
    
    BV_trans_init = BV0;
    
    i = 1;
    
    measurements_TN_FP = [];
    measurements = [];
    p = 1;
    
    if H~= 0
        while x1>-per_BV_critical*BV_trans_init && Vr>0 && Vr>0.1*p*(x1+BV0) 
        %while x1>-0.5 && Vr>0 && x1<(p)*BV_trans_init 
        
        if i<=trantime/st
            h = 0;
            u = 0;
            BV_trans_init = x1+BV0;
            p = 0.005;
            Vr_temp = Vr;
        else
            h = H;
            u = I;
            p = 1;
        end
        
        x1_dot_bar = [x1_dot_bar;-K*x1 + x2 + u - h];
        x2_dot_bar = [x2_dot_bar;K*u/(1+alpha_u) - K*h/(1+alpha_h)];
        
        x1 = x1 + st*(-K*x1 + x2 + u - h);
        x2 = x2 + st*(K*u/(1+alpha_u) - K*h/(1+alpha_h));
        Vr = Vr - st*h*Vr/(x1+BV0);
        
        prev_z = z;
        z = (HCT0-(Vr/(x1+BV0)))/HCT0;
        z_dot = [z_dot;(z-prev_z)/st];
        
        v = v + st*u;
        
        HCT_bar = [HCT_bar;Vr/(x1+BV0)];
        x1_bar = [x1_bar;x1];
        x2_bar = [x2_bar;x2];
        v_bar = [v_bar;v];
        
        i = i+1;
        
        end

    states.x1_true = x1_bar;
    states.x2_true = x2_bar;
    states.x1_dot_true = x1_dot_bar;
    states.x2_dot_true = x2_dot_bar;
    states.BV_true = x1_bar+BV0;
    measurements.v = v_bar;
    HCT_bar = add_noise(HCT_bar,noise,noise_type);
    
    

    if filter_measurements
        HCT_bar = lowpass(HCT_bar,0.5,1/st);
    end

   
    measurements.HCT = HCT_bar;
    y_true = (HCT0-HCT_bar)./HCT_bar;
    measurements.y_true = y_true;
    measurements.z_dot = z_dot;

    else
        while x1>-per_BV_critical*BV_trans_init && Vr>0 && Vr>0.1*p*(x1+BV0) 
        %while x1>-0.5 && Vr>0 && x1<(p)*BV_trans_init 
        if i<=trantime/st
            h = 0;
            u = 0;
            BV_trans_init = x1+BV0;
            p = 0.005;
            Vr_temp = Vr;
        else
            h = H;
            u = I;
            p = 1;
        end
        
        x1_dot_bar = [x1_dot_bar;-K*x1 + x2 + u - h];
        x2_dot_bar = [x2_dot_bar;K*u/(1+alpha_u) - K*h/(1+alpha_h)];
        
        x1 = x1 + st*(-K*x1 + x2 + u - h);
        x2 = x2 + st*(K*u/(1+alpha_u) - K*h/(1+alpha_h));
        Vr = Vr - st*h*Vr/(x1+BV0);
        
        prev_z = z;
        z = (HCT0-(Vr/(x1+BV0)))/HCT0;
        z_dot = [z_dot;(z-prev_z)/st];
        
        v = v + st*u;
        
        HCT_bar = [HCT_bar;Vr/(x1+BV0)];
        x1_bar = [x1_bar;x1];
        x2_bar = [x2_bar;x2];
        v_bar = [v_bar;v];
        i = i+1;
        
        end

        states.x1_bar = x1_bar;
        states.x2_bar = x2_bar;
        states.x1_dot_bar = x1_dot_bar;
        states.x2_dot_bar = x2_dot_bar;
        measurements_TN_FP.v = v_bar;
        HCT_bar = add_noise(HCT_bar,noise,noise_type);
        if filter_measurements
            HCT_bar = lowpass(HCT_bar,0.5,1/st);
        end
        measurements_TN_FP.HCT = HCT_bar;
        y_true = (HCT0-HCT_bar)./HCT_bar;
        measurements_TN_FP.y_true = y_true;
        measurements_TN_FP.z_dot = z_dot;
    
    end

    
       

end