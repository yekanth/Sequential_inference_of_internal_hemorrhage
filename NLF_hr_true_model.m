
function [y,HCT,x_T] = NLF_hr_true_model(VSUBJECTS,Hrate,Irate,st,ic)
    
    % Non-linear True or Actual Model
    Dt = st;
    x_0 = ic;
    
    x1 = x_0(1);
    x2 = x_0(2);
    
    
    JH = Hrate;
    JI = Irate;
    
    theta_s = VSUBJECTS;
    
    alpha_u = theta_s(:,1);
    alpha_h = theta_s(:,2);
    k      = theta_s(:,3);
    V0      = theta_s(:,11);
    H0      = theta_s(:,12);
    
    Vr = V0*H0;
    
    x_true = [];
    
    for i = 1:length(Irate)
        x1_next = x1 + Dt*(-k*x1 + x2 - JH(i) + JI(i));
        x2_next = x2 + Dt*((-k/(1+alpha_h))*JH(i) + (k/(1+alpha_u)*JI(i)));
        Vr_next = Vr - Dt*(JH(i)*(Vr/(x1+V0)));
        HCT(i) = Vr/(x1+V0);
        
        x_true = [x_true; [x1 x2 Vr]];
        
        x1 = x1_next;
        x2 = x2_next;
        Vr = Vr_next;
        
       
    end
    
    y = (H0 - HCT)./HCT;
    
    x_T = x_true;
    
end
