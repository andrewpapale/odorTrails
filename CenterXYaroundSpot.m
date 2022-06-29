% center x,y around spot.

nM = max(mouse);
nS = max(sess);


xbins = length(linspace(-1300,1300,round(1300*2/(11.2*5))));
ybins = length(linspace(-1100,1100,round(1100*2/(11.2*5))));
M = zeros(xbins,ybins);
Hxy = M;
Hnxy = M;
Hdt = M;
Hdnt = M;
Hv = M;
Hnv = M;
Hdphi = M;

for iM=1:nM
    for iS=1:nS
        k = mouse==iM & sess==iS;
        kT = mouseT==iM & sessT==iS;
        
        if sum(k)>0 && t2s(iM,iS)>1
            
            k = k & ~edge & t2s0==1;
            
            x0 = x(k);
            y0 = y(k);
            nx0 = nx(k);
            ny0 = ny(k);
            nV0 = nV(k);
            V0 = V(k);
            dphi0 = log10dphi(k);
            tempt2s0 = t2s0(k);
            tempspotfound = spotfound1(k);
            cx0 = nanmedian(xT1(kT));
            cy0 = nanmedian(yT1(kT));
            
            % center relative to spot
            x0 = x0-cx0;
            y0 = y0-cy0;
            nx0 = nx0-cx0;
            ny0 = ny0-cy0;
            
            % get aligned histograms
            H0 = histcn([x0,y0],xbins,ybins);
            
            % fix histogram weirdness with padding
            if size(H0,1)>size(Hxy,1)
                H0 = cat(1,H0(1:size(Hxy,1),:));
                %disp('truncating H0...');
            elseif size(H0,1)<size(Hxy,1)
                H0 = cat(1,H0,nan(1,size(Hxy(:,1))-size(H0(:,1))));
                %disp('padding H0...');
            end
            if size(H0,2)>size(Hxy,2)
                H0 = cat(2,H0(:,1:size(Hxy,2)));
                %disp('truncating H0...');
            elseif size(H0,2)<size(Hxy,2)
                H0 = cat(2,H0,nan(size(Hxy(:,2))-size(H0(:,2))));
            end           
            
            Hxy = squeeze(sum(cat(3,Hxy,H0),3));
           
            % get aligned histograms
            H1 = histcn([nx0,ny0],xbins,ybins);
            
            % fix histogram weirdness with padding
            if size(H1,1)>size(Hnxy,1)
                H1 = cat(1,H1(1:size(Hnxy,1),:));
                %disp('truncating H0...');
            elseif size(H1,1)<size(Hnxy,1)
                H1 = cat(1,H1,nan(1,size(Hnxy(:,1))-size(H1(:,1))));
                %disp('padding H0...');
            end
            if size(H1,2)>size(Hnxy,2)
                H1 = cat(2,H1(:,1:size(Hnxy,2)));
                %disp('truncating H0...');
            elseif size(H1,2)<size(Hnxy,2)
                H1 = cat(2,H1,nan(size(Hnxy(:,2))-size(H1(:,2))));
            end           
            
            Hnxy = squeeze(sum(cat(3,Hnxy,H1),3));
            
            % get aligned histograms
            H0 = histcn([nx0,ny0],xbins,ybins,'AccumData',dT(k),'fun',@nanmean);
            
            % fix histogram weirdness with padding
            if size(H0,1)>size(Hdt,1)
                H0 = cat(1,H0(1:size(Hdt,1),:));
                %disp('truncating H0...');
            elseif size(H0,1)<size(Hdt,1)
                H0 = cat(1,H0,nan(1,size(Hdt(:,1))-size(H0(:,1))));
                %disp('padding H0...');
            end
            if size(H0,2)>size(Hdt,2)
                H0 = cat(2,H0(:,1:size(Hdt,2)));
                %disp('truncating H0...');
            elseif size(H0,2)<size(Hdt,2)
                H0 = cat(2,H0,nan(size(Hdt(:,2))-size(H0(:,2))));
            end           
            
            Hdt = squeeze(nanmean(cat(3,Hdt,H0),3));
            
             % get aligned histograms
            H0 = histcn([nx0,ny0],xbins,ybins,'AccumData',dnT(k),'fun',@nanmean);
            H0(H1==0)=nan;
            % fix histogram weirdness with padding
            if size(H0,1)>size(Hdnt,1)
                H0 = cat(1,H0(1:size(Hdnt,1),:));
                %disp('truncating H0...');
            elseif size(H0,1)<size(Hdnt,1)
                H0 = cat(1,H0,nan(1,size(Hdnt(:,1))-size(H0(:,1))));
                %disp('padding H0...');
            end
            if size(H0,2)>size(Hdnt,2)
                H0 = cat(2,H0(:,1:size(Hdnt,2)));
                %disp('truncating H0...');
            elseif size(H0,2)<size(Hdnt,2)
                H0 = cat(2,H0,nan(size(Hdnt(:,2))-size(H0(:,2))));
            end           
            
            Hdnt = squeeze(nanmean(cat(3,Hdnt,H0),3));
            
            % get aligned histograms
            H0 = histcn([nx0,ny0],xbins,ybins,'AccumData',V(k),'fun',@nanmean);
            H0(H1==0)=nan;
            % fix histogram weirdness with padding
            if size(H0,1)>size(Hv,1)
                H0 = cat(1,H0(1:size(Hv,1),:));
                %disp('truncating H0...');
            elseif size(H0,1)<size(Hv,1)
                H0 = cat(1,H0,nan(1,size(Hv(:,1))-size(H0(:,1))));
                %disp('padding H0...');
            end
            if size(H0,2)>size(Hv,2)
                H0 = cat(2,H0(:,1:size(Hv,2)));
                %disp('truncating H0...');
            elseif size(H0,2)<size(Hv,2)
                H0 = cat(2,H0,nan(size(Hv(:,2))-size(H0(:,2))));
            end           
            
            Hv = squeeze(nanmean(cat(3,Hv,H0),3));
            
            % get aligned histograms
            H0 = histcn([nx0,ny0],xbins,ybins,'AccumData',nV(k),'fun',@nanmean);
            H0(H1==0)=nan;
            % fix histogram weirdness with padding
            if size(H0,1)>size(Hnv,1)
                H0 = cat(1,H0(1:size(Hnv,1),:));
                %disp('truncating H0...');
            elseif size(H0,1)<size(Hnv,1)
                H0 = cat(1,H0,nan(1,size(Hnv(:,1))-size(H0(:,1))));
                %disp('padding H0...');
            end
            if size(H0,2)>size(Hnv,2)
                H0 = cat(2,H0(:,1:size(Hnv,2)));
                %disp('truncating H0...');
            elseif size(H0,2)<size(Hnv,2)
                H0 = cat(2,H0,nan(size(Hnv(:,2))-size(H0(:,2))));
            end           
            
            Hnv = squeeze(nanmean(cat(3,Hnv,H0),3));

            % get aligned histograms
            H0 = histcn([nx0,ny0],xbins,ybins,'AccumData',log10dphi(k),'fun',@nanmean);
            H0(H1==0)=nan;
            % fix histogram weirdness with padding
            if size(H0,1)>size(Hdphi,1)
                H0 = cat(1,H0(1:size(Hdphi,1),:));
                %disp('truncating H0...');
            elseif size(H0,1)<size(Hnv,1)
                H0 = cat(1,H0,nan(1,size(Hdphi(:,1))-size(H0(:,1))));
                %disp('padding H0...');
            end
            if size(H0,2)>size(Hdphi,2)
                H0 = cat(2,H0(:,1:size(Hdphi,2)));
                %disp('truncating H0...');
            elseif size(H0,2)<size(Hdphi,2)
                H0 = cat(2,H0,nan(size(Hdphi(:,2))-size(H0(:,2))));
            end           
            
            Hdphi = squeeze(nanmean(cat(3,Hdphi,H0),3)); 
            
        end
        disp(iS);
    end
    disp(iM);
end
            