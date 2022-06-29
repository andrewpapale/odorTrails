function vid_subtractBackground_opto(videoPath)

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

Line =[536.39,12.56; 555.59, 37.82; 576.16,61.94];
result = computeline([Line(2,1),Line(2,2)],(Line(3,2)-Line(1,2))./(Line(3,1)-Line(1,1)), [520 580]);
nP = length(result);
x1 = result{1}(:,1);
y1 = result{1}(:,2);
x2 = result{nP}(:,1);
y2 = result{nP}(:,2);

nT = 100;
tbins = linspace(1,vr.Duration-0.1,nT);
background = [];
for iT=1:nT
    vr.currentTime = tbins(iT);
    I = vr.readFrame;
    I = double(I);
    [j,i]=meshgrid(1:size(I,1),1:size(I,2));
    %zeroOut = (i-576).^2+(j-0).^2 < 100.^2;
    Dline = ((y2-y1)*i-(x2-x1)*j+x2*y1-y2*x1)./sqrt((y2-y1).^2+(x2-x1).^2);
    zeroOut = Dline > 0;
    I(zeroOut')=0;
    I = uint8(I);
    background = cat(3,background,I(:,:,1));
    %background(background > 250)= nan;
end

medianbkg = nanmedian(background,3);


% get "mouse background"



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
net = denoisingNetwork('DnCNN'); 
Opto_LED_Signal = [];
frame = [];
iF = 1;



while ~done
    % Read next frame
    I = vr.readFrame();
    I = double(I);
    [j,i]=meshgrid(1:size(I,1),1:size(I,2));
    %zeroOut = (i-576).^2+(j-0).^2 < 100.^2;
    Dline = ((y2-y1)*i-(x2-x1)*j+x2*y1-y2*x1)./sqrt((y2-y1).^2+(x2-x1).^2);
    zeroOut = Dline > 0;
    I1 = I;
    I(zeroOut')=0;
    I1(~zeroOut')=nan;
    I = uint8(I);
    %I = adapthisteq(I(:,:,1)-medianbkg);
    %I = imsharpen(adapthisteq(imadjust(I,[0.01 1],[0 1],0.5),'NumTiles',[20 20],'ClipLimit',0.8,'Distribution','exponential','alpha',1),'Threshold',0.5);
    %I = adapthisteq(imsharpen(imadjust(I-medianbkg,[0.01,1],[0,1],0.9),'Threshold',0.8),'NumTiles',[80 80],'ClipLimit',0.1,'Distribution','uniform');
    %I = imadjust((I-medianbkg)-imgaussfilt(I-medianbkg,0.5,'FilterDomain','spatial'),[0.001,1],[0,1],0.33);
    %I = locallapfilt((I-medianbkg)-imgaussfilt(I-medianbkg,0.5,'FilterDomain','spatial'),0.05,0.3,'NumIntensityLevels',20);
    I = (I-0.9*medianbkg)-0.98*imgaussfilt(I-0.9*medianbkg,5,'FilterDomain','spatial');
    I = denoiseImage(I,net);
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
    
    Opto_LED_Signal = cat(1,Opto_LED_Signal,nanmean(I1(:)));
    frame = cat(1,frame,iF+1);
end
elapsed = toc(t0);
fprintf('Finished writing %d frames in %.2f mins.\n', framesWritten, elapsed/60)
close(V1);

savePath1 = strcat(saveStr{1},'-LED');
save(savePath1,'Opto_LED_Signal','frame');
