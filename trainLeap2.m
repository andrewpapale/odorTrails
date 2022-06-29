
D = dir('trainLEAPbox1.h5');
currdir = cd;
savePath = strcat(currdir,'\','trainLEAPbox.h5');

Dstr = D(1).name;
vidInfo = h5info(Dstr);
nX = vidInfo.Datasets.Dataspace.Size(1);
nY = vidInfo.Datasets.Dataspace.Size(2);
frames = vidInfo.Datasets.Dataspace.Size(4);

h5create(savePath,'/box',[nX,nY,1,inf],'ChunkSize',[nX,nY,1,1],'Deflate',1,'Datatype','uint8');

framesWritten = 0;
for iT=1:frames
    vid = h5read(Dstr,'/box1',[1 1 1 iT],[nX,nY,1,1]);
    h5write(savePath, '/box', vid, [1 1 1,framesWritten+1], [nX,nY,1,1]);
    framesWritten = framesWritten+1;
    disp(iT);
end