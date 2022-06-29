% load '/media/james/ExtraDrive1/Matlab/Code/AP Video Data/2018-03-19-AndyP.mat'
% load '/media/james/ExtraDrive1/Matlab/Code/AP Video Data/backup matlab/Spot Paper/Figure Scripts/model_trajectories.mat'

% gett2s0
% getspotfound1


%%
edge = ~(nearX>44.8 & nearY>44.8 & nearX<nanmax(nearX)-44.8 & nearY<nanmax(nearY)-44.8);
isstopped = zeros(size(x)); nM=5; nS=197; for iM=1:nM; for iS=1:nS; k=mouse==iM & sess==iS & ~edge; if sum(k)>0; [xd,~,~]=density_filter(nx(k),ny(k),20,0.2); isstopped(find(isnan(xd)))=1; end; end; disp(iM); end

%%
bin_edge = 0:140;
% Calculate midpoints of bins for plotting
mdpts_s = bin_edge(1:end-1)+(bin_edge(2)-bin_edge(1))/2;
%Compute area of each ring corresponding to each bin
ringArea = nan(1,length(bin_edge)-1);
for i = 2:length(bin_edge)-1
    ringArea(i) =  2*pi*bin_edge(i)^2 - 2*pi*bin_edge(i-1)^2 ;
end

%% Data

F = figure;

%Successful
k = t2s0 ==1 & spotfound1 ==1 & ~edge & ~isstopped;
% Pool ALL successful trial frames and bin by distance
N_s = histcounts(dnT(k),bin_edge);
% Compute frac occupancy
frac_occ = N_s./sum(N_s);
% Normalize by area in ring
frac_occ = frac_occ./ringArea;
% Plot log10(frac_occ)
plot(mdpts_s, log10(frac_occ), 'LineWidth' , 2, 'Color' , 'r')
hold on

%Unsuccessful
k = spotfound1 ==0 & ~edge & ~isstopped;
% Pool ALL unsuccessful trial frames and bin by distance
N_f = histcounts(dnT(k),bin_edge);
% Compute frac occupancy
frac_occ = N_f./sum(N_f);
% Normalize by area in ring
frac_occ = frac_occ./ringArea;
% Plot log10(frac_occ)
plot(mdpts_s, log10(frac_occ), 'LineWidth' , 2, 'Color' , [.5 .5 .5 ])
hold on

%% Model

%Cat all distances for successful trials
numRuns = length(mod_out);
dists_s = [];
for i = 1:numRuns
    if mod_out(i).success == 1
        dists_s = [dists_s ; mod_out(i).dists];
    end
end

% Pool ALL successful trial frames and bin by distance
N_m = histcounts(dists_s,bin_edge);
% Compute frac occupancy
frac_occ = N_m./sum(N_m);
% Normalize by area in ring
frac_occ = frac_occ./ringArea;
% Plot log10(frac_occ)
plot(mdpts_s, log10(frac_occ), 'LineWidth' , 2, 'Color' , 'k')
hold on



title('Fractional occupancy vs distance')
xlim([0 100])
xlabel('Distance from odor source (cm)')
ylabel('log_{10}( fractional occupancy / ring area) (cm^{-2})')
set(gca,'FontSize',18)