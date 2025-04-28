function BP = estimate_BP(SN,H,I,trantime,measurements,VSUBJECTS_NAT,nominal,noise,st,mt)
    
    % Function to estimate Blood Pressure using the mathematical model and
    % hemorrhage, resusciation inputs
    HCT_measured = measurements.HCT;
    theta = VSUBJECTS_NAT(SN,:);
    n_VP = length(VSUBJECTS_NAT(:,1));
    
    alltime = length(HCT_measured)*mt;
    
    theta_s = nominal;
    
    if length(HCT_measured)*mt>=trantime
        Irate = I.*[zeros(trantime/st,1); ones(ceil(alltime/st)-trantime/st,1)];
    else
        Irate = zeros(ceil(alltime/st),1);
    end

    if length(HCT_measured)*mt>=trantime
        Hrate = H.*[zeros(trantime/st,1); ones(ceil(alltime/st)-trantime/st,1)];
    else
        Hrate = zeros(ceil(alltime/st),1);
    end
    
    % Time
    ti=(length(Irate)-1)*st;
    T_end = ti; % End time [minutes]
    st = 0.01; % Sampling time [minutes]
    time_range = 0:st:T_end;
    k_end = length(time_range);
    
    % Inputs

    UI = Irate;
    UH = Hrate;

    U = [UI,UH];
    
     % Parameters
    alpha_Gain = theta(:,1);
    alpha_Loss = theta(:,2);
    Kp  = theta(:,3);
    Ka  = theta(:,4);
    Kv  = theta(:,5);
    a = theta(:,6);
    b = theta(:,7);
    Ktp = theta(:,8);
    beta_v = theta(:,9);
    Kco = theta(:,10);
    
    BV0 = theta(:,11);
    HCT0 = theta(:,12);
    CO0 = theta(:,13);
    BP0 = theta(:,14);
    
    CP0 = 6; % mmHg
    
    % Constants
    EQM = 1; % Remove inputs to find equilibrium
    AF  = 0.30; % Default fraction of blood volume in artery
    Va0 = AF*BV0;
    Vv0 = (1-AF)*BV0;

    R = (BP0-CP0)./CO0;
    BP_Target = BP0;
    
    Vr0 = HCT0.*(BV0);

    x0 = [Va0 Vv0 HCT0.*BV0 0 0 0];
    
    % Simulate model

    x = x0;
    X = zeros(length(x0), k_end);
    Y = zeros(5, k_end);
    
    
    for k = 1:k_end
        
        % Time
        t = time_range(k);

        % Save state
        X(:,k) = x;
    
        % Perform transition
        [ x, y ] = BP_transition_model( x, theta, U(k,:), st );

        % Save output
        BP(k) = y(3);
    
    end
    
    


end