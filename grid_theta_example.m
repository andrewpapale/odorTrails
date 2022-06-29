% These are x,y coordinates for the vectors
xs = linspace(-70,70,101);
ys = linspace(-70,70,101);


% These arrays store the average and standard deviation for each vector
% coordinate
theta_means = nan(length(xs),length(ys));
theta_stds = nan(length(xs),length(ys));

tic
for i = 1:length(xs)
    parfor j = 1:length(ys)
        % This loop retuns a list of heading coordinates th from
        % a model, and computes their mean and standard deviation
        
        [~,~,~,th] = klinotaxis_Beta_particle(xs(i),ys(j),2*pi*rand,'vel',0, 'Time',1000,'Capture',0 ,'AB',1);
        
        theta_means(i,j) = circ_mean(th);
        theta_stds(i,j) = circ_std(th);
        disp(num2str([i,j]))
        
    end
end
toc


% This section converts coordinates to list for plotting via quiver and scatter 
x = [];
y = [];
thm = [];
ths = [];

for i = 1:length(xs)
    for j = 1:length(ys)
        x = [x ; xs(i)];
        y = [y ; ys(j)];
        thm = [thm ; theta_means(i,j)];
        ths = [ths ; theta_stds(i,j)];
        
    end
end


% This plots the vectors at coords x,y.  Vector length direction is the
% mean direction of the heading th.  Vector length is inversely
% proportional to heading standard deviation.

quiver(x,y,cos(thm)./ths, sin(thm)./ths, .5, 'Color','k', 'LineWidth',2)
