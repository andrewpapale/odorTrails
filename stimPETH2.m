% stimPETH2
nM = max(mouse);
nS = max(sess);
dF = 100;
dt = 2;
tbins = 0.1;
%PETHt = nan(size(-dt:tbins:dt));
%PETHt = nan(length(-dt:tbins:dt),25);
PETHt = [];
PETHx = nan(size(-100:1:100));
Nt = [];
for iM=1:nM
    for iS=1:nS
        k = mouse==iM & sess==iS;
        kL = mouseL==iM & sessL==iS;
        if sum(k)>0
            %nV0 = zC(k);
            nV0 = nV(k);
            %nV0 = orient(k);
            %nV0 = log10dphi(k);
            t0 = t(k);
            frame0 = frame(k);
            dnT0 = dnT(k);
            k1 = edge(k) | nV(k) > 100 | bV(k) > 100 & t2s2(k)==0;
            nV0(k1)=nan;
            timeOut0 = timOut(kL);
            fraOut0 = find(stimF(k)==1);
            nP = sum(kL);
            zeroframe = [];
            for iP=1:nP
                zeroframe = interp1(t0,frame0,timeOut0(iP),'nearest');
                t00 = unique(interp1(t0,t0,timeOut0(iP)-dt:0.02:timeOut0(iP)+dt,'nearest'));
                f11 = unique(interp1(t0,frame0,t00,'nearest'));
                t00(isnan(f11))=[];
                f11(isnan(f11))=[];
                nV11 = nV0(f11);
                dnT11 = dnT0(f11);
                center = find(f11==zeroframe);
                %disp(sum(~isnan(nV11)));
%                 Pt0 = nan(size(PETHt));
                 H = histcn(t00',(timeOut0(iP)-dt:tbins:timeOut0(iP)+dt),'AccumData',nV11,'fun',@nanmean);
                  H(end+1)=nan;
                  N = histcn(t00',(timeOut0(iP)-dt:tbins:timeOut0(iP)+dt));
                  N(end+1)=nan;
                  H(N<1)=nan;

                % orientation 2D
%               if sum(~isnan(nV11))>10
                    %H = histcn([t00',dnT11],(timeOut0(iP)-dt:tbins:timeOut0(iP)+dt),linspace(2,100,25),'AccumData',nV11,'fun',@nanmean);
                    %H = histcn([t00',dnT11],(timeOut0(iP)-dt:tbins:timeOut0(iP)+dt),linspace(2,100,25),'AccumData',nV11*pi/180,'fun',@circ_mean);
                    %H = H*180/pi;
                    %PETHt = cat(1,PETHt,H');
                    %H(end+1,:)=nan;
%                     if size(H,2)<25
%                         H(:,end+1)=nan;
%                     else
%                     end
%                H(H==0)=nan;
%                else
%                   H = nan(41,25);
                    %H = nan(61,25);
 %              end
                PETHt = cat(2,PETHt,H);
                %PETHt = cat(3,PETHt,H);
                %Nt = cat(1,Nt,N');
                
%                 fpre = find((f11-zeroframe)<0);
%                 lenbefore = length(fpre);
%                 Px0 = nan(201,1);
%                 Px0(center:length(f11))=nV11(center:length(f11));
%                 left = sum(isnan(Px0(1:center-1)));
%                 Px0(left-lenbefore+1:center-1)=nV11(fpre);
%                 
%                 PETHx = cat(1,PETHx,Px0');
%                 
                
                
            end
        end
        disp(iS);
    end
    disp(iM);
end

%PETHt(:,nansum(Nt)<10)=nan;
            