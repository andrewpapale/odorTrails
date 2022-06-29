function v_est = pi_diff(y, Ts, kp, ki, postSmoothing)
% 2017-07-03 AndyP added smoothing

error_acc = 0;
y_hat = 0;

y0 = y(~isnan(y));

nS = max(size(y0));
for k = 1 : nS
    error = y0(k) - y_hat;
    error_acc = nansum(cat(1,error_acc,error * Ts));
    v_est(k) = kp * error + ki * error_acc;
    y_hat = nansum(cat(1,y_hat,v_est(k)*Ts));   
end

v_est0 = nan(size(y));
v_est0(~isnan(y)) = v_est;
v_est = v_est0;

if postSmoothing
	nS = ceil(postSmoothing/Ts);
	v_est = conv2(v_est,ones(nS)/nS,'same');
end