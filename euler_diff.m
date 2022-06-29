function v_est = euler_diff(y, Ts, postSmoothing)
% 2017-07-03 AndyP added smoothing



y0 = y(~isnan(y));
v_est = zeros(1, length(y0));

nS = max(size(y0));
for k = 1 : nS - 1
    v_est(k + 1) = (y0(k + 1) - y0(k)) / Ts;
end

v_est0 = nan(size(y));
v_est0(~isnan(y)) = v_est;
v_est = v_est0;

if postSmoothing
	nS = ceil(postSmoothing/Ts);
	v_est = conv2(v_est,ones(nS)/nS,'same');
end