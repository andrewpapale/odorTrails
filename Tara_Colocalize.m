% Tara_Colocalize
% 2019-11-05-AndyP
% written to colocalize PV and VVA

fprintf('open merged image \n');
% open image
[FileName,~,~] = uigetfile;
I = imread(FileName);
fprintf('open PV file \n');
uiopen;
areaPV = areas1;
posPV = cell_pos1;
radiusPV = sqrt(areaPV./pi);
fprintf('open VVA file \n');
uiopen;
radiusVVA = radius;
posVVA = cell_pos1;

figure(1); clf;
imagesc(I);
hold on;
%viscircles(posPV,radiusPV,'color','r');
viscircles(posVVA,radiusVVA,'color','g');
plot(posPV(:,1),posPV(:,2),'r.','markersize',10);
axis ij;
set(gca,'fontsize',18);
axis off;

% % determine if radius of PV cell overlaps with any radius of VVA cell
% th = 0:pi/50:2*pi;
% nPV = size(posPV,1); % number of PV cells counted
% nVVA = size(posVVA,1);
% nCoLo = nan(nPV,nVVA);
% for iPV = 1:nPV
%     for iVVA=1:nVVA
%         [xout,yout] = circcirc(posPV(iPV,1),posPV(iPV,2),radiusPV(iPV),posVVA(iVVA,1),posVVA(iVVA,2),radiusVVA(iVVA));
%  

% determine if PV center of mass is within VVA circle
nPV = size(posPV,1);
nVVA = size(posVVA,1);
ixCoLo = nan(nPV,nVVA);
nCoLo = 0;
for iVVA=1:nVVA
    ix = find((posPV(:,1)-posVVA(iVVA,1)).^2+(posPV(:,2)-posVVA(iVVA,2)).^2 <= radiusVVA(iVVA).^2);
    if ~isempty(ix)
        ixCoLo(ix,iVVA)=1;
        nCoLo = nCoLo+1;
        plot(posPV(ix,1),posPV(ix,2),'ms','markersize',10);
        if length(ix)>1
            warning('multiple cells colocalized with VVA cell %d',iVVA);
        end
    end
end

dateStr = strsplit(datestr(now),' ');
tempStr = strsplit(FileName,'_');
saveStr = strcat(tempStr{1},'_','colocalization','_',dateStr{1},'.mat');
save(saveStr,'nCoLo','ixCoLo');
