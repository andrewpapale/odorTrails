% make LEAP video
% 2019-02-14 AndyP

nF = 5000;
modelPath = 'C:\Users\papalea\Documents\Data\New Spot with Automatic Feeder\models\190214_161508-n=41\final_model.h5';
vid = h5read('356074-2018-07-26-1.h5','/box',[1,1,1,1],[480,576,1,nF]);


preds = predict_box(vid, modelPath);

bx = squeeze(preds.positions_pred(1,1,:));
by = squeeze(preds.positions_pred(1,2,:));
nx = squeeze(preds.positions_pred(2,1,:));
ny = squeeze(preds.positions_pred(2,2,:));
tbx = squeeze(preds.positions_pred(3,1,:));
tby = squeeze(preds.positions_pred(3,2,:));
ttx = squeeze(preds.positions_pred(4,1,:));
tty = squeeze(preds.positions_pred(4,2,:));

%%

V1 = VideoWriter('2019-02-14-testVideo-annotated','Uncompressed AVI');
set(V1,'FrameRate',round(1./nanmedian(diff(time0))));
open(V1);
for iF=500:nF
    F = figure(1);
    imagesc(squeeze(vid(:,:,1,iF))); 
    colormap gray;
    axis xy;
    hold on;
    plot(nx0(iF),ny0(iF),'g.','markersize',10);
    plot(nx(iF),ny(iF),'co','markersize',10);
    plot(x0(iF),y0(iF),'rx','markersize',10);
    plot(bx(iF),by(iF),'mo','markersize',10);
    viscircles([cx cy],pixpercm*radiusAroundSpot);
    title(strcat(mat2str(iF),'/',mat2str(nF)));
    frame = getframe(gcf);
    writeVideo(V1,frame);
    clf;  
end
close(V1);