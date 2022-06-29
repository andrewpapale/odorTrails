newV = VideoWriter('downsample-2.avi'); 
iF=1;
newV.FrameRate = V.FrameRate/3; 
open(newV);
nF = V.FrameRate*V.Duration;
for iF=1:nF
    if mod(iF,3)==0
        disp(iF); 
        currFrame = read(V,iF); 
        writeVideo(newV,currFrame); 
    end 
end
close(newV);