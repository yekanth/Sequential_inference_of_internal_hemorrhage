function [NDT,DT] = calculate_NDT_last(I,tran_time,window,prob_threshold,T_critical,st)
    
    %Do not consider anything before Transient Time 
    I = round(I(1:end),0);
    P = zeros(length(I),1);
    %Calculate the probability of Hemorrhage
    for i = tran_time/st+1:length(I)
        % if i<=(tran_time/st)+(window/st)
        %     P(i) = length(find(I(tran_time/st:i)==1))/i;
        % else
        if i>ceil(window/st)+tran_time/st
            P(i) = (length(find(I(i-ceil(window/st):i)==1)))*st/window;
        else
            P(i) = (length(find(I(tran_time/st:i)==1)))/(i-(tran_time/st)+1);
        end
        %end
    end
    %Calculate when the probability of hemorrhage crosses threshold and
    %remains above threshold
    above_threshold = 0;
    crossindices = [];
    if ~isempty(P)
        for k = tran_time/st+1:length(P)
            if P(k)>prob_threshold && ~above_threshold
                crossindices = [crossindices k];
                above_threshold = 1;
            elseif P(k)<prob_threshold
                above_threshold = 0;
            end
        end
    end



    if ~isempty(P)&&above_threshold
        PP = P;
        new_window = window;
        % for j = new_window/st+1:length(P)
        %     PP(j) = length(find(P(j-floor(new_window/st):j)>prob_threshold))>prob_threshold*new_window/st;
        % end
            
        T_det = min(find(PP>=prob_threshold))*st;
    else
        T_det = [];
    end
        
    
    NDT = (T_det-tran_time)/(T_critical-tran_time);
    DT = T_det;
end