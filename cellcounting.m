cc=bwconncomp(cellbw2);
centers = regionprops(cc,'Centroid');
volume = regionprops3(cc,'Volume');
radius = regionprops3(cc,'PrincipalAxisLength');

iz = nan(cc.NumObjects,1);
ix = nan(cc.NumObjects,1);
iy = nan(cc.NumObjects,1);
for iC=1:cc.NumObjects
    ix(iC)=centers(iC).Centroid(1);
    iy(iC)=centers(iC).Centroid(2);
    iz(iC)=round(centers(iC).Centroid(3));
end

% vi = nan(cc.NumObjects,1);
% for iC=1:cc.NumObjects
%     vi(iC)=volume(iC);
% end

ri = table2array(radius);
A = pi*ri(:,1)/2.*ri(:,2)/2.*0.321; % um ^2
vi = table2array(volume);


%%
for iZ=1:100
    subplot(1,2,1);
    imagesc(I2(:,:,iZ));
    subplot(1,2,2);
    imagesc(wat(:,:,iZ));
    hold on;
    if any(iz==iZ)
        plot(ix(iz==iZ),iy(iz==iZ),'rx');
    end
    pause;
end


