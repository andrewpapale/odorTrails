

npars = 10000;
iters = 100;


rng('shuffle');
kbin    = 0   + 300*rand(npars,1);      % [0 300]

rng('shuffle');
k       = 0.01 + 0.99*rand(npars,1);    % [0.01 1]

rng('shuffle');
noiselevel = rand(npars,1);             % [0 1]

rng('shuffle');
filt = 10*rand(npars,1);                % [0 10]








pars = [kbin k  noiselevel filt];

successes = nan(npars,1);




for i = 1:length(pars)
    tic
    successlocal = nan(iters,1);
    parfor j = 1:iters
        [~,~,~,~,~,~,success] = concTrackFun_smooth_comp_ou_active_nav_paper_version_corr(100,0,1.5,...
            'Time', 30,...
            'kbin',kbin(i),...
            'Vmax', 25,...
            'L', 5, ...
            'dnares', 0.18 ,...
            'phi_min',0.2 ,...
            'phi_max',0.3,...
            'k', k(i),...
            'noiselevel', noiselevel(i),...
            'filt',filt(i));
        
        successlocal(j) = success;
    end
    successes(i) = mean(successlocal);
    
    
    if mod(i,100)==0
        disp([num2str(i) ' of ' num2str(length(pars)) '.  Fraction complete: ' num2str(i/length(pars)) '. Success rate: ' num2str(successes(i))])
        toc
    end
    
end

save('2019-07-29-screen_reduced.mat','successes','pars');

