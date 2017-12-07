% get orientation from spot

orient = nan(size(x));
borient = nan(size(x));
dO = nan(size(x));
% IF = nan(size(x));
% IA = nan(size(x));
% IP = nan(size(x));

welch = [];
welch2 = [];
mouse2 = [];
S = [];
freqs = linspace(0.1,10,100);
params.tapers = [3 5];
params.Fs = 50;
%params.fpass = linspace(0.01,10,250);
mouse3 = [];
S1 = [];
nM = max(mouse);
nS = max(sess);
miss = 0;
spotfound = [];
Oezt = [];
Azdphi = [];
Sdnt = [];

for iM=1:nM
    for iS=1:nS
        k = mouse==iM & sess==iS;
        
        if sum(k)>0
            
            kT = mouseT==iM & sessT==iS;
            
            xT = nanmedian(xT1(kT));
            yT = nanmedian(yT1(kT));
            
            x0 = x(k);
            y0 = y(k);
            nx0 = nx(k);
            ny0 = ny(k);
            nearX0 = nearX(k);
            nearY0 = nearY(k);
            
            %// Index array for factor
            x1 = 1:numel(x0);
            %// Indices of NaNs
            t2 = find(~isnan(x0));
            
            %// Replace NaNs with the closest non-NaNs
            x0 = interp1(x1(t2),x0(t2),x1,'linear')';
            y0 = interp1(x1(t2),y0(t2),x1,'linear')';
            
            
            
            %// Index array for factor
            nx1 = 1:numel(nx0);
            %// Indices of NaNs
            t2 = find(~isnan(nx0));
            
            %// Replace NaNs with the closest non-NaNs
            nx0 = interp1(nx1(t2),nx0(t2),nx1,'linear')';
            ny0 = interp1(nx1(t2),ny0(t2),nx1,'linear')';
            
            A = atan2(yT-y0,xT-x0);
            B = atan2(ny0-y0,nx0-x0);
            orient(k)=angdiff(A,B)*180/pi;
            
            C = atan2(y0,x0);
            D = atan2(ny0,nx0);
            borient(k)=angdiff(C,D)*180/pi;
            
            dO(k) = foaw_diff(orient(k),1/50,50,0.5,0.1);
            
            t2s0 = t2s(iM,iS);
            
            if t2s0>0
                t2s0 = round(t2s0*50);
                spotfound = cat(1,spotfound,1);
                lastspotfound = t2s0;
            elseif t2s0==-20
                t2s0 = lastspotfound;
                spotfound = cat(1,spotfound,0);                
            elseif t2s0==-40
                t2s0 = lastspotfound;
                spotfound = cat(1,spotfound,-1);
            end
            if t2s0 > 50
                welch0 =pwelch(orient(1:t2s0),[],[],freqs,50);
                welch = cat(1,welch,welch0);
                mouse2 = cat(1,mouse2,repmat(iM,[1,length(freqs)]));
                
                
                welch1 = pwelch(dO(:),[],[],freqs,50);
            else
                spotfound = spotfound(1:end-1);
            end
            
            %[S0,times,freqs1]=mtspecgramc(borient(k),[1 0.05],params);
            
            %subplot(1,3,1);
            %pcolor(times,freqs1,S0'); shading flat; colorbar;
            %pause;
            
%             zS = [];
%             for iF=1:size(S0,2)
%                 zS(:,iF)=nanzscore(S0(:,iF));
%             end
            
            dnT0 = dnT(k);
            
            edge = find(~(nearX0>50 & nearY0>50 & nearX0<1280-50 & nearY0<1024-50));
            dnT0(edge) = nan;
            
            
            %tinsec = linspace(0,sum(k)/50,sum(k));
            %frame = interp1(tinsec,1:sum(k),times,'nearest');
            
            %Sdnt0 = dnT0(frame);
            %Sdnt = cat(1,Sdnt,Sdnt0);
            %S0(find(isnan(Sdnt0)),:)=nan; %#ok<FNDSB>
            %S = cat(1,S,S0);
            
            %subplot(1,3,2);
            %hist(Sdnt0,30);
            
             %sortSbydnT
%             
             %subplot(1,3,3);
             %pcolor(linspace(0,90,nB),freqs1,Ssort'); shading flat;
             %pause(0.5);
            
        end
    end
    disp(iM);
end