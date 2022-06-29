%cx = nan(size(x)); cy = nan(size(x)); for iM=1:nM; for iS=1:nS; k = mouse==iM & sess==iS; if sum(k)>0; cx(k) = x(k)-cx0(iM,iS); cy(k)=y(k)-cy0(iM,iS); end; end; end;



% center on stim
dt = 3;
cx = [];
cy = [];
torientc = [];
thetac = [];
orientc = [];
dphic = [];
for iM=1:nM
    for iS=1:nS
        k = mouse==iM & sess==iS;
        if sum(k)>0
            k1 = ~edge(k) & bV(k) < 100 & nV(k) < 100 & t2s2(k)==1;
            kL = mouseL==iM & sessL==iS;
            nx0 = nx(k);
            ny0 = ny(k);
            nV0 = nV(k);
            torient0 = torient(k);
            orient0 = orient(k);
            dphi0 = log10dphi(k);
            theta0 = theta(k);
            nx0(k1)=nan;
            ny0(k1)=nan;
            t0 = t(k);
            frame0 = frame(k);
            timOut0 = timOut(kL);
            fraOut0 = fraOut(kL);
            for iT=1:sum(kL)
                zeroframe = interp1(t0,frame0,timOut0(iT),'nearest');
                t00 = unique(interp1(t0,t0,timOut0(iT):0.015:timOut0(iT)+dt,'nearest'));
                f11 = unique(interp1(t0,frame0,timOut0(iT):0.015:timOut0(iT)+dt,'nearest'));
                f11(isnan(f11))=[];
                center = find(f11==zeroframe);
                
                xc0 = nx0(f11)-nx0(fraOut0(iT));
                yc0 = ny0(f11)-ny0(fraOut0(iT));
                
                f12 = unique(interp1(t0,frame0,timOut0(iT)-1:0.015:timOut0(iT),'nearest'));
                f12(isnan(f12))=[];
                
                alpha = nanmean(orient0(f12));
                
                
%                 M = [xc0, yc0]';
%                 R = [cosd(alpha), -sind(alpha); sind(alpha), cosd(alpha)];
%                 M1 = R*M;
%                 nxcr = M1(1,:);
%                 nycr = M1(2,:);

                cx = cat(1,cx,xc0);
                cy = cat(1,cy,yc0);
                torientc = cat(1,torientc,torient0(f11));
                thetac = cat(1,thetac,theta0(f11));
                orientc = cat(1,orientc,orient0(f11));
                dphic = cat(1,dphic,dphi0(f11));
                

                
            end
        end
    end
end











xbins = linspace(-300,300,100);
ybins = linspace(-300,300,80);

iC = 1;
X = [];
Y = [];
for iT=1:length(xbins)-1
    for iQ=1:length(ybins)-1
        X(iC) = xbins(iT);
        Y(iC)= ybins(iQ);
        iC = iC+1;
    end
end


%k = ~edge & ~isstopped & t2s2==1 & nV < 100 & V < 100;

H = histcn([cx,cy],xbins,ybins,'AccumData',(orientc*pi/180),'fun',@circ_mean);
dH = histcn([cx,cy],xbins,ybins,'AccumData',(orientc*pi/180),'fun',@circ_std);
dPhi = histcn([cx,cy],xbins,ybins,'AccumData',dphic,'fun',@nanmean);
dH(dH==0)=nan;
N = histcn([cx,cy],xbins,ybins);
H(N==0)=nan;
H = H';
H = H(:);
dPhi = dPhi';
dPhi = dPhi(:);

F = figure(1); clf;
k0 = abs(sin(H(:))./dH(:)) < pi & abs(cos(H(:))./dH(:)) < pi;
q = quiver(X(k0)'./5.174,Y(k0)'./5.174,cos(H(k0))./dH(k0),sin(H(k0))./dH(k0),2,'LineWidth',0.8,'MaxHeadSize',15,'AlignVertexCenters','on');
%// Compute the magnitude of the vectors
mags = dPhi(k0);
%// Get the current colormap
currentColormap = colormap(gca);

%// Now determine the color to make each arrow using a colormap
%[~, ~, ind] = histcounts(mags, size(currentColormap, 1));
%set(gca,'clim',[0 10]);
set(gca,'clim',[-2 0]);
clims = num2cell(get(gca, 'clim'));
[~, ~, ind] = histcounts(mags, linspace(clims{:}, size(currentColormap, 1)));

%// Now map this to a colormap to get RGB
cmap = uint8(ind2rgb(ind(:), currentColormap) * 255);
cmap(:,:,4) = 255;
cmap = permute(repmat(cmap, [1 3 1]), [2 1 3]);

%// We repeat each color 3 times (using 1:3 below) because each arrow has 3 vertices
set(q.Head, ...
    'ColorBinding', 'interpolated', ...
    'ColorData', reshape(cmap(1:3,:,:), [], 4).');   %'

%// We repeat each color 2 times (using 1:2 below) because each tail has 2 vertices
set(q.Tail, ...
    'ColorBinding', 'interpolated', ...
    'ColorData', reshape(cmap(1:2,:,:), [], 4).');

hold on;
viscircles([0,0],1.5);
axis off;
colorbar;
set(gca,'fontsize',18);
