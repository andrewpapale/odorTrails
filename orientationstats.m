notedge = x > 44.8 & x < (1280-44.8) & y > 44.8 & y < (1024-44.8);
nB = 30;
bdnT = quantile(dnT(notedge),nB);

p02 = nan(nB,1);
k02 = nan(nB,1);
K02 = nan(nB,1);
p12 = nan(nB,1);
k12 = nan(nB,1);
K12 = nan(nB,1);

h02 = nan(nB,1);
h12 = nan(nB,1);

% p0 = nan(nB,1);
% p1 = nan(nB,1);
% p2 = nan(nB,1);
% k0 = nan(nB,1);
% k1 = nan(nB,1);
% k2 = nan(nB,1);

for iT=1:nB
    if iT==nB
        keep = dnT>bdnT(iT);
    else
        keep = dnT>bdnT(iT) & dnT<=bdnT(iT+1);
    end
    
    k00 = keep & notedge & nV>0.1 & t2s0==0 & hitmiss1==0; % spot missed
    k11 = keep & notedge & nV>0.1 & t2s0==0 & hitmiss1==1; % after spot found
    k22 = keep & notedge & nV>0.1 & t2s0==1 & hitmiss1==1; % before spot found
    
    %
    %     H0 = hist(orient(k00),linspace(-180,180,360));
    %     H1 = hist(orient(k11),linspace(-180,180,360));
    %     H2 = hist(orient(k22),linspace(-180,180,360));
    %
    %     H0 = H0./nansum(H0);
    %     H1 = H1./nansum(H1);
    %     H2 = H2./nansum(H2);
    
    %     figure(1); clf;
    %     plot(linspace(-180,180,360),H0); hold on;
    %     plot(linspace(-180,180,360),H1);
    %     plot(linspace(-180,180,360),H2);
    
    %    pause;
%     if sum(k00)>0 && sum(k11)>0 && sum(k22)>0
%         figure(1); clf;
%         cdfplot(orient(k00));
%         hold on;
%         cdfplot(orient(k11));
%         cdfplot(orient(k22));
%         legend('missed','after found','before found');
%         title(mat2str(iT));
%         pause;
%     end
    
    if sum(k00)>100 && sum(k22)>100
        %[p0(iT),k0(iT)]=circ_raotest(orient(k00)*pi/180);
        [p02(iT),k02(iT),K02(iT)] = circ_kuipertest(orient(k00)*pi/180,orient(k22)*pi/180,360,0);
        %h02(iT)= p02(iT) < 0.05/nB;
        
        %[h02(iT),p02(iT),k02(iT)]=kstest2(orient(k00),orient(k22),'Alpha',0.001/nB);
    end
    
    if sum(k11)>100 && sum(k22)>100
        %[p1(iT),k1(iT)]=circ_raotest(orient(k11)*pi/180);
        [p12(iT),k12(iT),K12(iT)] = circ_kuipertest(orient(k11)*pi/180,orient(k22)*pi/180,360,0);
        %h12(iT)= p12(iT) < 0.05/nB;
        
        %[h12(iT),p12(iT),k12(iT)]=kstest2(orient(k11),orient(k22),'Alpha',0.001/nB);
    end
    
    %     if sum(k22)>10
    %         [p2(iT),k2(iT)]=circ_raotest(orient(k22)*pi/180);
    %     end
end

% figure(1); clf;
% plot(h02); hold on;
% plot(h12);
%
% figure(2); clf;
% plot(log10(p02)); hold on;
% plot(log10(p12));


