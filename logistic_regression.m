

close all; clc

k = t2s>-40 & abs(meandphi)<1000;

Y = spotfound(k);
var_1 =  ALorKP0(k);

% figure
% plot(Y, 'x')


succ_ind = Y == 1;
fall_ind = Y == 0;

figure(1); hold on
plot(var_1(succ_ind(:)), ones(sum(succ_ind(:))), 'gx')
plot(var_1(fall_ind(:)), zeros(sum(fall_ind(:))), 'rx')



%%

% for iM=1:nM
%     for iS=1:nS
% for iM=1
%     for iS=1:30
% %         k = mouse==iM & sess==iS;
%         k = mouse==iM & sess==iS & spotfound(iM, iS) == 0;
% %         kT = mouseT==iM & sessT==iS;
%         if sum(k)>0
% %             xT = nanmedian(xT1(kT));
% %             yT = nanmedian(yT1(kT));
% %             x0 = nx(k);
% %             y0 = ny(k);
%             dist_tmp = dnT(k);
%             figure
%             plot(dist_tmp)
% 
% %             R1 = 
% %             dist_tmp = dnT(k);
% %             find(dist_tmp > )
%         
%         end
%             
%     end
% end



%% logistic regression
% 
% x = [2100 2300 2500 2700 2900 3100 ...
%      3300 3500 3700 3900 4100 4300]';
% n = [48 42 31 34 31 21 23 23 21 16 17 21]';
% y = [1 2 0 3 8 17 14 8 19 15 17 21]';
% yy = (y./n)>0.5;
% 
% b = glmfit(x, yy, 'binomial', 'link', 'probit');
% 
% yfit = glmval(b, x,'probit');
% figure
% plot(x, yy,'o',x, yfit,'-','LineWidth',2)














