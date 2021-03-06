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
params.tapers = [3 9];
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
            
%             %// Index array for factor
%             x1 = 1:numel(x0);
%             %// Indices of NaNs
%             t2 = find(~isnan(x0));
%             
%             %// Replace NaNs with the closest non-NaNs
%             x0 = interp1(x1(t2),x0(t2),x1,'spline')';
%             y0 = interp1(x1(t2),y0(t2),x1,'spline')';
%             
%             
%             
%             %// Index array for factor
%             nx1 = 1:numel(nx0);
%             %// Indices of NaNs
%             t2 = find(~isnan(nx0));
%             
%             %// Replace NaNs with the closest non-NaNs
%             nx0 = interp1(nx1(t2),nx0(t2),nx1,'spline')';
%             ny0 = interp1(nx1(t2),ny0(t2),nx1,'spline')';
            
            A = atan2(yT-y0,xT-x0);
            B = atan2(ny0-y0,nx0-x0);
            orient(k)=angdiff(A,B)*180/pi;
            
            C = atan2(y0,x0);
            D = atan2(ny0,nx0);
            borient(k)=angdiff(C,D)*180/pi;
            
%             dO(k) = foaw_diff(orient(k),1/50,50,0.5,0.1);
%             
%             
%             welch0 =pwelch(orient(k),[],[],freqs,50);
%             welch = cat(1,welch,welch0);
%             mouse2 = cat(1,mouse2,repmat(iM,[1,length(freqs)]));
%             
%             
%             welch1 = pwelch(dO(k),[],[],freqs,50);
%             welch2 = cat(1,welch,welch1);
%             
%             % need to align to entering 15cm times
%             [S0,times,freqs1]=mtspecgramc(borient(k),[2 0.1],params);
%             
%             dnT0 = dnT(k);
%             ezt = find(dnT0<=15,1,'first');
%             Oezt = cat(1,Oezt,orient(ezt));
%             
%             if ~isempty(ezt)
%                 
%                 dt = nanmedian(diff(times));
%                 t1 = (ezt-400)/50;
%                 t2 = (ezt+400)/50;
%                 
%                 times1 = find(times >= t1 & times <= t2);
%                 timesfull = t1:round(dt,2,'significant'):t2;
%                 lfull = length(timesfull);
%                 edge = find(~(nearX0>50 & nearY0>50 & nearX0<1280-50 & nearY0<1024-50));
%                 edge(edge<times1(1) | edge > times1(end))=[];
%                 goodtimes = find(timesfull>0 & timesfull<max(times))+times1(1)-1;
%                 for iE=1:length(edge)
%                     goodtimes(goodtimes==edge(iE))=[];
%                 end
%                 goodtimes(goodtimes>max(times1))=[];
%                 
%                 if ~isempty(goodtimes)
%                     times2 = interp1(times1,times1,goodtimes);
%                     
%                     S1 = S0(times2,:);
%                     if ezt < 400
%                         pad = nan(lfull-size(S1,1),size(S0,2));
%                         S1 = cat(1,pad,S1);
%                     elseif ezt > size(S0,1)-400
%                         pad = nan(lfull-size(S1,1),size(S0,2));
%                         S1 = cat(1,S1,pad);
%                     else
%                         goodtimes = linspace(goodtimes(1),goodtimes(end),lfull);
%                         times2 = round(interp1(times1,times1,goodtimes),1,'significant');
%                         S1 = S0(times2,:);
%                     end
%                     
%                     
%                     if ~isempty(S1)
%                         try
%                             S = cat(3,S,S1);
%                             mouse3 = cat(1,mouse3,iM);
%                             if t2s(iM,iS)>0
%                                 spotfound = cat(1,spotfound,1);
%                             elseif t2s(iM,iS)==-20
%                                 spotfound = cat(1,spotfound,0);
%                             elseif t2s(iM,iS)==-40
%                                 spotfound = cat(1,spotfound,-1);
%                             end
%                         catch
%                             keyboard;
%                         end
%                     end
%                 else
%                     miss = miss+1;
%                     disp('all times along edge');
%                 end
%                 
%                 
%                 % align zdphi1
%                 edge = find(~(nearX0>50 & nearY0>50 & nearX0<1280-50 & nearY0<1024-50));
%                 
%                     S1 = S0(times2,:);
%                     if ezt < 400
%                         pad = nan(lfull-size(S1,1),size(S0,2));
%                         S1 = cat(1,pad,S1);
%                     elseif ezt > size(S0,1)-400
%                         pad = nan(lfull-size(S1,1),size(S0,2));
%                         S1 = cat(1,S1,pad);
%                     else
%                         goodtimes = linspace(goodtimes(1),goodtimes(end),lfull);
%                         times2 = round(interp1(times1,times1,goodtimes),1,'significant');
%                         S1 = S0(times2,:);
%                     end
%                     
%                 
%                 
%                 
%             else
%                 miss = miss+1;
%                 disp('did not cross 15cm');
%             end
            %[IF(k), IA(k), IP(k)]=InstFreq(orient(k));
        end
    end
    disp(iM);
end