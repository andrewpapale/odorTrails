function [v_est,winused] = foaw_diff(y, Ts, m, d, postSmoothing)
% 2017-07-03 AndyP added smoothing


slope = 0;                                  % estimate


v_est0 = nan(size(y));
winused = nan(size(y));
y0 = y(~isnan(y));
nS = max(size(y0));

if ~isempty(y0)
    for k = 2 : nS
        window_len = 0;
        can_increase_window = true;
        
        while 1
            window_len = window_len + 1;
            
            if (window_len > m || k - window_len == 0)
                window_len = window_len - 1;
                break;
            end
            
            % slope of the line of: y(k) = slope * k * Ts + c
            % this line is passing through y_k to y_k-i
            slope_ = slope;
            slope = (y0(k) - y0(k - window_len)) / (window_len * Ts);
            
            if (window_len > 1)
                c = y0(k) - slope * k * Ts;
                
                % Check every point from k to k-j
                for j = 1 : window_len - 1
                    delta = y0(k - j) - (c + slope * (k - j) * Ts);
                    if (abs(delta) > 2*d)
                        can_increase_window = false;
                        window_len = window_len - 1;
                        slope = slope_;
                        break;
                    end %% end if
                end %% end for
                
            end %% end if
            
            if can_increase_window == false
                break;
            end
        end
        winused(k) = window_len;
        v_est(k) = slope;
    end
    
    if postSmoothing
        nS = ceil(postSmoothing/Ts);
        v_est = conv2(v_est,ones(nS)/nS,'same');
    end
    
    v_est0(~isnan(y))=v_est;
end

v_est = v_est0;



