function data = getData_getTrailGUI(handles,nbrh)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

minThr = str2num(get(handles.edit1,'String')); %#ok<ST2NM>
maxThr = str2num(get(handles.edit2,'String')); %#ok<ST2NM>
sigma = round(str2num(get(handles.edit3,'String'))); %#ok<ST2NM>
Thr = str2num(get(handles.edit5,'String'));  %#ok<ST2NM>
Noise = str2num(get(handles.edit6,'String')); %#ok<ST2NM>
edgeVal = round(str2num(get(handles.edit4,'String'))); %#ok<ST2NM>
edgeDetection = get(handles.checkbox1, 'Value');
medianFilter = get(handles.checkbox2,'Value');
Filling = get(handles.checkbox3,'Value');
ThresholdData = get(handles.checkbox4,'Value');
if medianFilter && ~edgeDetection && ~ThresholdData
    data = medfilt2(handles.current_data,[nbrh nbrh]);
elseif medianFilter && edgeDetection && ~ThresholdData
    data = medfilt2(edge(handles.current_data,'Canny',[minThr maxThr],sigma),[nbrh nbrh]);
elseif ~medianFilter && edgeDetection && ~ThresholdData
    data = edge(handles.current_data,'Canny',[minThr maxThr],sigma);
elseif ~medianFilter && ~edgeDetection && ~ThresholdData
    data = handles.current_data;
elseif ThresholdData
    [y,x]=find(handles.current_data>Thr*nanmean(handles.current_data(:)));
    data = zeros(size(handles.current_data));
    for iX=1:length(x)
        data(y(iX),x(iX))=1;
    end
    
    if medianFilter
        data = medfilt2(data,[nbrh,nbrh]);
    else
    end
end

data(1:edgeVal,:)=0;
data(end-edgeVal:end,:)=0;
data(:,1:edgeVal)=0;
data(:,end-edgeVal:end)=0;

CC = bwconncomp(data);
numPixels = cellfun(@numel,CC.PixelIdxList);
idx = find(numPixels<Noise);
for i=1:length(idx)
    data(CC.PixelIdxList{idx(i)})=0;
end

if Filling && edgeDetection
    data = ConenctClosestPoints(data);
    data = FillingAlgorithm(data);
    data = bwmorph(data,'thin',round(sigma));
elseif Filling && ThresholdData
    data = ConenctClosestPoints(data);
else
end

end

