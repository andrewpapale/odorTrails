function Wrap_LEAP

modelPath = 'C:\Users\papalea\Documents\Data\New Spot with Automatic Feeder\models\190221_135044-n=61\final_model.h5';
%modelPath = 'C:\Users\papalea\Documents\Data\Hartung Data\models\190416_095642-n=51';


Ddone = dir('*.h5');
chunkSize = 3000;

for iD=1:length(Ddone)
    %try
    DdoneStr = Ddone(iD).name;
    saveStr = strsplit(DdoneStr,'.');
    vid = [];
    if ~exist(strcat(saveStr{1},'.mat'))
        vidInfo = h5info(DdoneStr);
        nX = vidInfo.Datasets.Dataspace.Size(1);
        nY = vidInfo.Datasets.Dataspace.Size(2);
        frames = vidInfo.Datasets.Dataspace.Size(4);
        
        bx = [];
        by = [];
        nx = [];
        ny = [];
        tbx = [];
        tby = [];
        ttx = [];
        tty = [];
        confb = [];
        confn = [];
        conftb = [];
        conftt = [];
        
        frameBins = ceil(linspace(0,frames,frames/chunkSize+1));
        nB = length(frameBins);
        for iT=1:nB-1
            if iT<nB-1
                vid = h5read(DdoneStr,'/box',[1 1 1 frameBins(iT)+1],[nX,nY,1,frameBins(iT+1)-frameBins(iT)]);
            else
                vid = h5read(DdoneStr,'/box',[1 1 1 frameBins(iT)+1],[nX,nY,1,frames-(frameBins(iT)+1)]);
            end  
            preds = predict_box(vid, modelPath);
            bx = cat(1,bx,squeeze(preds.positions_pred(1,1,:))); %#ok<*NASGU>
            by = cat(1,by,squeeze(preds.positions_pred(1,2,:)));
            nx = cat(1,nx,squeeze(preds.positions_pred(2,1,:)));
            ny = cat(1,ny,squeeze(preds.positions_pred(2,2,:)));
            tbx = cat(1,tbx,squeeze(preds.positions_pred(3,1,:)));
            tby = cat(1,tby,squeeze(preds.positions_pred(3,2,:)));
            ttx = cat(1,ttx,squeeze(preds.positions_pred(4,1,:)));
            tty = cat(1,tty,squeeze(preds.positions_pred(4,2,:)));
            confb = cat(1,confb,squeeze(preds.conf_pred(1,:))');
            confn = cat(1,confn,squeeze(preds.conf_pred(2,:))');
            conftb = cat(1,conftb,squeeze(preds.conf_pred(3,:))');
            conftt = cat(1,conftt,squeeze(preds.conf_pred(4,:))');
        end
        
        
        save(saveStr{1},'bx','by','nx','ny','confb','confn','ttx','tty','tbx','tby','conftb','conftt');
        disp(saveStr);
        %catch
        %    warning('could not process video %s',DdoneStr);
        %end
    end
end