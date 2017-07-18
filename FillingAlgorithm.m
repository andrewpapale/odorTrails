function L = FillingAlgorithm(data)

% obtained on 2017-06-13 from URL:
% https://www.mathworks.com/matlabcentral/answers/23467-filling-the-objects-found-by-canny-edge
cc = bwconncomp(data);
% As you can see, some connected components are not *closed*:
%figure, imagesc(imfill(labelmatrix(cc),'holes')), drawnow
% So, let's try closing them by iterative dilations
BWblank = false(cc.ImageSize);
stats = regionprops(cc,'ConvexImage','EulerNumber');
for i = find([stats.EulerNumber]>0)
    distIm = bwdist(~stats(i).ConvexImage);
    maxClose = ceil(max(distIm(:)));
    BWslice = BWblank;
    BWslice(cc.PixelIdxList{i}) = true;
    if isinf(maxClose), continue; end;
    for dilSz = 2:maxClose
        BWnew = imdilate(BWslice,ones(dilSz));
        statsNew = regionprops(BWnew,'EulerNumber');
        if statsNew.EulerNumber<=0
            BWnew = imerode(imfill(BWnew,'holes'),ones(dilSz));
            cc.PixelIdxList{i} = find(BWnew);
        end
    end
end
%figure, imagesc(imfill(labelmatrix(cc),'holes')), drawnow
% That got almost all of them. Some are left over where the dilation itself
% filled everything so the euler number stayed at 1. Let's just replace
% those with their convex hull
stats = regionprops(cc,'ConvexImage','EulerNumber','BoundingBox');
for i = find([stats.EulerNumber]>0)
    maxClose = ceil(max(distIm(:)));
    BWslice = BWblank;
    BWslice(cc.PixelIdxList{i}) = true;
    distIm = bwdist(~BWslice);
    if ~any(distIm(:)>1)
        BWnew = BWslice;
        bb = ceil(stats(i).BoundingBox);
        BWnew((1:bb(4))+bb(2)-1,(1:bb(3))+bb(1)-1) = stats(i).ConvexImage;
        cc.PixelIdxList{i} = find(BWnew);
    end
end
L = imfill(labelmatrix(cc),'holes');
%figure, imagesc(L)
% Now we know that any blobs surrounded by other blobs are actually holes
indsOfHoles = find(arrayfun(@(i)mode(double(L(bwmorph(L==i,'dilate',1)&~(L==i)))),1:cc.NumObjects));
L(ismember(L,indsOfHoles)) = 0;
% Hint: get(hObject,'Value') returns toggle state of checkbox3

L = bwareaopen(L,100); % remove spots with fewer than 100 pixels
end

