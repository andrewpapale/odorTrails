function vid_subtractBackground(videoPath)

%% Parameters
% Path to input file
%videoPath = 'C:\Users\papalea\Documents\Data\New Spot with Automatic Feeder\336712-2018-06-27-2.avi';
currDir = cd;
videoPath = strcat(currDir,'\',videoPath);
% Path to output file
% savePath = '..\..\leap\data\examples\072212_163153.clip.h5';
%savePath = 'C:\Users\papalea\Documents\Data\New Spot with Automatic Feeder\336712-2018-06-27-2.h5';
saveStr = strsplit(videoPath,'.');
savePath = strcat(saveStr{1},'-bkg','.avi');
% Frames to convert at a time (lower this if your memory is limited)
chunkSize = 1000;

% Convert frames to single channel grayscale images (instead of 3 channel RGB)
grayscale = true;

%% Initialize
% Open video for reading
vr = VideoReader(videoPath);

% Check size from first frame
I0 = vr.readFrame();
if grayscale
    %I0 = imadjust(rgb2gray(I0)); 
    %I0 = imadjust(I0); 
end
frameSize = size(I0);
if numel(frameSize) == 2; frameSize = [frameSize 1]; end

nT = 100;
tbins = linspace(1,vr.Duration-0.1,nT);
background = [];
for iT=1:nT
    vr.currentTime = tbins(iT);
    I = vr.readFrame;
    background = cat(3,background,I(:,:,1));
    %background(background > 250)= nan;
end

medianbkg = nanmedian(background,3);

% Reset VideoReader
delete(vr);
vr = VideoReader(videoPath);

% Check if file already exists
if exist(savePath,'file') > 0
    warning(['Overwriting existing bkg subtracted file: ' savePath])
    delete(savePath)
end

%% Save
buffer = cell(chunkSize,1);
done = false;
framesRead = 0;
framesWritten = 0;
t0 = tic;
V1 = VideoWriter(savePath,'Grayscale AVI');
V1.FrameRate = 35;
open(V1);
while ~done
    % Read next frame
    I = vr.readFrame();
    %I = adapthisteq(I(:,:,1)-medianbkg);
    %I = imsharpen(adapthisteq(imadjust(I,[0.01 1],[0 1],0.5),'NumTiles',[20 20],'ClipLimit',0.8,'Distribution','exponential','alpha',1),'Threshold',0.5);
    %I = adapthisteq(imsharpen(imadjust(I-medianbkg,[0.01,1],[0,1],0.9),'Threshold',0.8),'NumTiles',[80 80],'ClipLimit',0.1,'Distribution','uniform');
    I = imadjust((I-medianbkg)-imgaussfilt(I-medianbkg,1,'FilterDomain','spatial'),[0.001,1],[0,1],0.33);
    %I = imadjust(I,[0.25 0.95])-medianbkg;
    % Check if there are any frames left
    F = figure(1); clf;
    imagesc(I);
    set(gcf,'Position',[0 0 576 480]);
    axis off;
    done = ~vr.hasFrame();
    % Increment frames read counter and add to the write buffer
    framesRead = framesRead + 1;
    disp(framesRead);
    F1 = getframe(F);
    F1.cdata = F1.cdata(:,:,1);
    writeVideo(V1,F1);
end
elapsed = toc(t0);
fprintf('Finished writing %d frames in %.2f mins.\n', framesWritten, elapsed/60)
close(V1);
