% 2017-07-13 AndyP
% This script centers spots and computes the average measure 'measure' 

nM = length(unique(mouse));
nS = max(sess);
window = 200; % pixels
window = window./5.174;
nB = 20;

H = [];
N = [];
for iM=1:nM
    for iS=1:nS
        k0 = mouse==iM & sess==iS;
        kT = mouseT==iM & sessT==iS;
        
        nx0 = nx(k0);
        ny0 = ny(k0);
        dnT0 = dnT(k0);
        znidphi0 = log10dphi(k0);
        %znC0 = znC(k0);
        
        xT0 = cx0(iM,iS);
        yT0 = cy0(iM,iS);
        
        kW = nx0 >= xT0-window & nx0 <= xT0+window & ny0 >= yT0-window & ny0 <= yT0+window;
         
%         figure(2); clf;
%         plot(nx0,ny0,'k.')
%         hold on
%         scatter(nx0,ny0,[],dnT(k0)); colorbar; caxis([0 30]);
%         scatter(nx0(kW),ny0(kW),[],dnT0(kW),'filled'); 
%         plot(xT1(kT),yT1(kT),'r.');
%       
%        pause;
        
        if sum(kW)>0
            
            
            
            
            H0 = histcn([nx0(kW),ny0(kW)],linspace(xT0-window,xT0+window,nB),linspace(yT0-window,yT0+window,nB),'AccumData',znidphi0(kW),'fun',@nanmean);
            H0 = padarray(H0,[nB-size(H0,1),nB-size(H0,2)],nan,'post');
            N0 = histcn([nx0(kW),ny0(kW)],linspace(xT0-window,xT0+window,nB),linspace(yT0-window,yT0+window,nB));
            N0 = padarray(N0,[nB-size(N0,1),nB-size(N0,2)],nan,'post');
            H0(N0==0)=nan;
            H = cat(3,H,H0);
            N = cat(3,N,N0);
        end
    end
end


