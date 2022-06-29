function [e, eimage] = detectEdgesInFrame(filename, time)
% function e = detectEdges(I)
%
% Detects the edges within an movie frame

vid_struct = mmread(filename, [], [time time+1]);
gf = rgb2gray(vid_struct.frames(1).cdata);
ei = edge(gf, 'canny');
props = {'Area', 'PixelIdxList'};
e = regionprops(ei, props);
ea = [e.Area];
long = ea >= 50;
e = e(long);
ea = ea(long);
[~, order] = sort(ea);
e = e(order);
e = e(end-10:end);

ei = zeros(size(ei));
for ii=1:length(e)
    ei(e(ii).PixelIdxList) = 1;
end
eimage = (ei==1);