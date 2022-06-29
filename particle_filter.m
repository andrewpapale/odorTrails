function[x_est_out,y_est_out] = particle_filter(nx0, ny0,tx0,ty0)
%Student Dave's particle filter tutorial

%Adapted from Dan Simon Optimal state estimation book and Gordon, Salmond,
%and Smith paper.


%% initialize the variables


start = find(~isnan(nx0),1,'first');
x = nx0(start); % initial actual state
y = ny0(start);
x_N = 2; % Noise covariance in the system (i.e. process noise in the state update, here, we'll use a gaussian.)
x_R = 10; % Noise covariance in the measurement (i.e. the Quail creates complex illusions in its trail!)
T = length(nx0); % duration the chase (i.e. number of iterations).
N = 40; % The number of particles the system generates. The larger this is, the better your approximation, but the more computation you need.

%initilize our initial, prior particle distribution as a gaussian around
%the true initial value

V = 1; %define the variance of the initial esimate
x_P = nan(N,1); % define the vector of particles
y_P = nan(N,1);
% make the randomly generated particles from the initial prior gaussian distribution
for i = 1:N
    x_P(i) = x + sqrt(V) * randn;
    y_P(i) = y + sqrt(V) * randn;
end

%{
%show the distribution the particles around this initial value of x.
figure(1)
clf
subplot(121)
plot(1,x_P,'.k','markersize',5)
xlabel('time step')
ylabel('flight position')
subplot(122)
hist(x_P,100)
xlabel('flight position')
ylabel('count')
pause
%}

%the functions used by the Quail are:
% x = 0.5*x + 25*x/(1 + x^2) + 8*cos(1.2*(t-1)) + PROCESS NOISE --> sqrt(x_N)*randn
% z = x^2/20 + MEASUREMENT NOISE -->  sqrt(x_R)*randn;

%generate the observations from the randomly selected particles, based upon
%the given function
%z_out = [x^2 / 20 + sqrt(x_R) * randn];  %the actual output vector for measurement values.
zx_out = abs(tx0(start)-nx0(start));
zy_out = abs(ty0(start)-ny0(start));
x_out = [nx0(1:start)];  %the actual output vector for measurement values.
x_est = [nx0(1:start)]; % time by time output of the particle filters estimate
x_est_out = [x_est]; % the vector of particle filter estimates.
y_out = [ny0(1:start)];
y_est = [ny0(1:start)];
y_est_out = [y_est];

dtail = cat(1,nan,sqrt(diff(tx0).^2+diff(ty0).^2));
dnose = cat(1,nan,sqrt(diff(nx0).^2+diff(ny0).^2));

ok = dtail < 20 & dnose < 20;
okidx = find(ok==1);
tidx = 1:T;

for t = start+1:T
    %from the previou time step, update the flight position, and observed
    %position (i.e. update the Quails position with the non linear function
    %and update from this position what the chasing ninja's see confounded
    %by the Quails illusions.
    x = nx0(t-1);
    y = ny0(t-1);
    zx = abs(nx0(t)-tx0(t))./(nx0(t)+tx0(t))+sqrt(x_R)*randn;
    zy = abs(ny0(t)-ty0(t))./(ny0(t)+ty0(t))+sqrt(x_R)*randn;
    
    %k2 = find(okidx > tidx(t-2) & okidx < tidx(k1),2,'last');

    %Here, we do the particle filter
    for i = 1:N
        %given the prior set of particle (i.e. randomly generated locations
        %the quail might be), run each of these particles through the state
        %update model to make a new set of transitioned particles.
        
        
        
        x_P_update(i) = sqrt(x_P(i)*nx0(t))+sqrt(x_R)*randn;
        zx_update(i) = x_P_update(i)*(abs(nx0(t)-tx0(t))./(nx0(t)+tx0(t)));
        y_P_update(i) = sqrt(y_P(i)*ny0(t))+sqrt(x_R)*randn;
        zy_update(i) = y_P_update(i)*(abs(ny0(t)-ty0(t))./(ny0(t)+ty0(t)));
        %with these new updated particle locations, update the observations
        %for each of these particles.
        
        %Generate the weights for each of these particles.
        %The weights are based upon the probability of the given
        %observation for a particle, GIVEN the actual observation.
        %That is, if we observe a location z, and we know our observation error is
        %guassian with variance x_R, then the probability of seeing a given
        %z centered at that actual measurement is (from the equation of a
        %gaussian)
        P_wx(i) = (1/sqrt(2*pi*x_R)) * exp(-(zx - zx_update(i))^2/(2*x_R));
        P_wy(i) = (1/sqrt(2*pi*x_R)) * exp(-(zy - zy_update(i))^2/(2*x_R));
    end
    
    % Normalize to form a probability distribution (i.e. sum to 1).
    P_wx = P_wx./nansum(P_wx);
    P_wy = P_wy./nansum(P_wy);
    %{
    figure(1)
    clf
    subplot(121)
    plot(P_w,z_update,'.k','markersize',5)
    hold on
    plot(0,z,'.r','markersize',50)
    xlabel('weight magnitude')
    ylabel('observed values (z update)')
    subplot(122)
    plot(P_w,x_P_update,'.k','markersize',5)
    hold on
    plot(0,x,'.r','markersize',50)
    xlabel('weight magnitude')
    ylabel('updated particle positions (x P update)')
    pause
    
    
    %plot the before and after
    figure(1)
    clf
    subplot(131)
    plot(0,x_P_update,'.k','markersize',5)
    title('raw estimates')
    xlabel('fixed time point')
    ylabel('estimated particles for flight position')
    subplot(132)
    plot(P_w,x_P_update,'.k','markersize',5)
    hold on
    plot(0,x,'.r','markersize',40)
    xlabel('weight magnitude')
    title('weighted estimates')
    %}
    %% Resampling: From this new distribution, now we randomly sample from it to generate our new estimate particles
    
    %what this code specifically does is randomly, uniformally, sample from
    %the cummulative distribution of the probability distribution
    %generated by the weighted vector P_w.  If you sample randomly over
    %this distribution, you will select values based upon there statistical
    %probability, and thus, on average, pick values with the higher weights
    %(i.e. high probability of being correct given the observation z).
    %store this new value to the new estimate which will go back into the
    %next iteration
    if ~all(isnan(P_wx))
        for i = 1 : N
            x_P(i) = x_P_update(find(rand <= nancumsum(P_wx),1));
        end
    end
    if ~all(isnan(P_wy))
        for i = 1 : N
            y_P(i) = y_P_update(find(rand <= nancumsum(P_wy),1));
        end
    end
    %The final estimate is some metric of these final resampling, such as
    %the mean value or variance
    x_est = nanmean(x_P);
    y_est = nanmean(y_P);
    %{
    %the after
    subplot(133)
    plot(0,x_P_update,'.k','markersize',5)
    hold on
    plot(0,x_P,'.r','markersize',5)
    plot(0,x_est,'.g','markersize',40)
    xlabel('fixed time point')
    title('weight based resampling')
    pause
    %}
    % Save data in arrays for later plotting
    x_out = cat(1,x_out,x);
    zx_out = cat(1,zx_out,zx);
    x_est_out = cat(1,x_est_out,x_est);
    y_out = cat(1,y_out,y);
    zy_out = cat(1,zy_out,zy);
    y_est_out = cat(1,y_est_out,y_est);
    
end


t = 1:T;
figure(1);
clf
subplot(1,2,1);
plot(t, x_out, '.-b', t, x_est_out, '-.r','linewidth',3);
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('time step'); ylabel('Quail flight position');
legend('True flight position', 'Particle filter estimate');
subplot(1,2,2);
plot(t, y_out, '.-b', t, y_est_out, '-.r','linewidth',3);
set(gca,'FontSize',12); set(gcf,'Color','White');
xlabel('time step'); ylabel('Quail flight position');
legend('True flight position', 'Particle filter estimate');


