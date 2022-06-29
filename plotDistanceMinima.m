%% Preliminaries

% load '/media/james/ExtraDrive1/Matlab/Code/AP Video Data/2018-03-19-AndyP.mat'
% load '/media/james/ExtraDrive1/Matlab/Code/AP Video Data/backup matlab/Spot Paper/Figure Scripts/model_trajectories.mat'
% gett2s0

promThres = 4;
%% Data

peaksM = [];
for m = 1:5
    for s = 1:197
        k = mouse == m & sess == s & t2s0 == 1 & V >0.1;
        if sum(k) > 0
            [pks,~,~,p] = findpeaks(-dnT(k) , 'MinPeakProminence',promThres );	% pks = peak magnitude, p = prominence measure, negative dnT to find minima
            peaksM = [peaksM;-pks];
        end
    end
end

% peaksM(find(isnan(peaksM)))=[];

%% Model

peaks = [];
for i = 1:2000
    if mod_out(i).success==1
        [pks,~,~,p] = findpeaks(-mod_out(i).dists, 'MinPeakProminence',promThres );
        peaks = [peaks;-pks];
    end
end

% peaks(find(isnan(peaks)))=[];


%% Plot


Hmod = histogram(peaks,linspace(0,140,100),'Normalization','probability')
hold on
Hdat = histogram(peaksM,linspace(0,140,100),'Normalization','probability')

axis tight
set(gca,'FontSize',18)
title('Distribution of nose distance-to-spot minima preceding capture')
xlabel('Nose distance-from-spot minima (cm)')
legend('Model','Data')