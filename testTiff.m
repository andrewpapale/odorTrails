% read tiff stack



for k = 1:100
 I(:,:,k) = imread('420723_2_2_021519_ROI2_2X_2umStep.tif', k);
 imagesc(I(:,:,k));
 pause;
end

