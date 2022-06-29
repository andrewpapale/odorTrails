% create training set

fpvid = 300;
D = dir('*.h5');
nD = length(D);

currdir = cd;

savePath = strcat(currdir,'\','trainLEAPbox1.h5');

load('2018-11-15-Dataset.mat','x','y','edge','mouse','sess');

iMtestSet = [ones(9,1);ones(10,1)+1;ones(8,1)+2];
iStestSet = [1,2,3,4,7,8,9,10,11,...
    1:10,...
    2:5,6:9]';

Dstr = D(1).name;
vidInfo = h5info(Dstr);
nX = vidInfo.Datasets.Dataspace.Size(1);
nY = vidInfo.Datasets.Dataspace.Size(2);

h5create(savePath,'/box1',[nX,nY,1,inf],'ChunkSize',[nX,nY,1,1],'Deflate',1,'Datatype','uint8');

framesWritten = 0;

for iD=1:nD
    Dstr = D(iD).name;
    vidInfo = h5info(Dstr);
    frames = vidInfo.Datasets.Dataspace.Size(4);

    iM = iMtestSet(iD);
    iS = iStestSet(iD);
    
    k = find(mouse==iM & sess==iS & ~edge & ~isnan(x))-find(mouse==iM & sess==iS,1,'first')+1;
    F = figure(1); clf;
    plot(x(k),y(k),'.');
    pause(0.1);
    %rframe = randperm(k(length(k)));
    %rframe = rframe(1:fpvid);
    
    for iR=1:min(fpvid,length(k))
        vid = h5read(Dstr,'/box',[1 1 1 k(iR)],[nX,nY,1,1]);
        h5write(savePath, '/box1', vid, [1 1 1,framesWritten+1], [nX,nY,1,1]);
        framesWritten = framesWritten+1;
    end
    fprintf('%d/%d \n',iD,nD);
end