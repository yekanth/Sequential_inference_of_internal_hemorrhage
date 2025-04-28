function [bar_max,min_sim_time,max_sim_time] = calculate_bar_max(vp_actual_measurements,vp_bar_states,n_VP,st,mt)

                %Find Maximum Actual Simulation Time Across all VPs
                max_sim_time = 0;
                min_sim_time = 1e8;
                for i=1:n_VP
                    field = sprintf('subject_%d', i);
                    sim_time = length(vp_actual_measurements.(field).HCT)*mt;
                    if sim_time>max_sim_time
                        max_sim_time = sim_time;
                    end
                    if sim_time<min_sim_time
                        min_sim_time = sim_time;
                    end
                end
                
                %Combine all the bar states into one array and zero pad
                for j = 1:max_sim_time/st
                    for i = 1:n_VP
                        field = sprintf('subject_%d', i);
                        this_VP_time = length(vp_actual_measurements.(field).HCT)*mt;
                            
                            if j<=floor(this_VP_time/st)
                                x1_bar_comp(i) = vp_bar_states.(field).x1_bar(j);
                                x1_dot_bar_comp(i) = vp_bar_states.(field).x1_dot_bar(j);
                                x2_bar_comp(i) = vp_bar_states.(field).x2_bar(j);
                                x2_dot_bar_comp(i) = vp_bar_states.(field).x2_dot_bar(j);
                            else
                                x1_bar_comp(i) = vp_bar_states.(field).x1_bar(floor(this_VP_time/st));
                                x1_dot_bar_comp(i) = vp_bar_states.(field).x1_dot_bar(floor(this_VP_time/st));
                                x2_bar_comp(i) = vp_bar_states.(field).x2_bar(floor(this_VP_time/st));
                                x2_dot_bar_comp(i) = vp_bar_states.(field).x2_dot_bar(floor(this_VP_time/st));
                            end
                    end
                    x1_bar_max(j) = select_percentile(x1_bar_comp,0.1,'descend');
                    x2_bar_max(j) = select_percentile(x2_bar_comp,0.1,'descend');
                    x1_dot_bar_max(j) = select_percentile(x1_dot_bar_comp,0.1,'descend');
                    x2_dot_bar_max(j) = select_percentile(x2_dot_bar_comp,0.1,'descend');
                end
                
                bar_max.x1_bar = x1_bar_max;
                bar_max.x2_bar = x2_bar_max;
                bar_max.x1_dot_bar = x1_dot_bar_max;
                bar_max.x2_dot_bar = x2_dot_bar_max;
                
                %Sanity Check Plot for bar states
                % figure(1)
                % subplot(2,1,1)
                % for i=1:n_VP
                %     field = sprintf('subject_%d', i);
                %     plot(vp_bar_states.(field).x1_bar(1:floor(min_sim_time/st)))
                %     hold on;
                % end
                % hold on;
                % plot(x1_bar_max(1:floor(min_sim_time/st)),'o')
                % 
                % subplot(2,1,2)
                % for i=1:n_VP
                %     field = sprintf('subject_%d', i);
                %     plot(vp_bar_states.(field).x2_bar(1:floor(min_sim_time/st)))
                %     hold on;
                % end
                % hold on;
                % plot(x2_bar_max(1:floor(min_sim_time/st)),'o')

end