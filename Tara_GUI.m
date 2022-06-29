function varargout = Tara_GUI(varargin)
%TARA_GUI MATLAB code file for Tara_GUI.fig
%      TARA_GUI, by itself, creates a new TARA_GUI or raises the existing
%      singleton*.
%
%      H = TARA_GUI returns the handle to a new TARA_GUI or the handle to
%      the existing singleton*.
%
%      TARA_GUI('Property','Value',...) creates a new TARA_GUI using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to Tara_GUI_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      TARA_GUI('CALLBACK') and TARA_GUI('CALLBACK',hObject,...) call the
%      local function named CALLBACK in TARA_GUI.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Tara_GUI

% Last Modified by GUIDE v2.5 05-Nov-2019 11:57:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Tara_GUI_OpeningFcn, ...
    'gui_OutputFcn',  @Tara_GUI_OutputFcn, ...
    'gui_LayoutFcn',  [], ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Tara_GUI is made visible.
function Tara_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for Tara_GUI
handles.output = hObject;
handles.Rect = [];
handles.Rect2 = [];
handles.Rect0 = [];
handles.Radius = [];
handles.tempRect1 = [];
[handles.FileName,~,~] = uigetfile;
handles.I = imread(handles.FileName);
handles.I0 = handles.I(:,:,1);
% handles.thrVals = [0.05 0.1 0.15];
% handles.sizeLims = [50 140];
handles.thrVal1 = 2.5;
handles.minSizeVal = 40;
handles.SensVal = 0.9;
handles.EdgeThrVal = 0.1;
handles.areas1 = [];
handles.eccentricity1 = [];
handles.solidity1 = [];
handles.areas2 = [];
handles.eccentricity2 = [];
handles.solidity2 = [];
handles.cell_pos2 = [];
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Tara_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Tara_GUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Exit.
function Exit_Callback(hObject, eventdata, handles)
% hObject    handle to Exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close all force

% --- Executes on button press in drawRect.
function drawRect_Callback(hObject, eventdata, handles)
% hObject    handle to drawRect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Rect0 = getrect;
handles.Rect = rectangle('Position',handles.Rect0);
set(handles.Rect,'EdgeColor','r','LineWidth',3);
handles.Area = handles.Rect.Position(3)*handles.Rect.Position(4);
set(handles.rectArea,'String',mat2str(handles.Area));
guidata(hObject,handles);

function rectArea_Callback(hObject, eventdata, handles)
% hObject    handle to rectArea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rectArea as text
%        str2double(get(hObject,'String')) returns contents of rectArea as a double


% --- Executes during object creation, after setting all properties.
function rectArea_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rectArea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SaveData.
function SaveData_Callback(hObject, eventdata, handles)
% hObject    handle to SaveData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Rect = handles.Rect;
Rect2 = handles.Rect2;
areas1 = handles.areas1;
cell_pos1 = handles.cell_pos1;
eccentricity1 = handles.eccentricity1;
solidity1 = handles.solidity1;
threshold = handles.thrVal1;
minSize = handles.minSizeVal;
radius = handles.Radius;
SensVal = handles.SensVal;
EdgeThrVal = handles.EdgeThrVal;


if handles.countAll.Value
    areas2 = [];
    cell_pos2 = [];
    eccentricity2 = [];
    solidity2 = [];
else
    areas2 = handles.areas2;
    cell_pos2 = handles.cell_pos2;
    eccentricity2 = handles.eccentricity2;
    solidity2 = handles.solidity2;
end
dateStr = strsplit(datestr(now),' ');
tempStr = strsplit(handles.FileName,'_');
saveStr = strcat(tempStr{1},'_','cellcounts','_',dateStr{1},'.mat');
save(saveStr,'SensVal','EdgeThrVal','radius','Rect','Rect2','areas1','cell_pos1','areas2','cell_pos2','eccentricity1','eccentricity2','solidity1','solidity2','threshold','minSize');

% --- Executes on button press in CopyRect.
function CopyRect_Callback(hObject, eventdata, handles)
% hObject    handle to CopyRect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(handles.Rect)
    handles.tempRect = imrect(gca,handles.Rect.Position);
end
while ~handles.LockRect2.Value
    guidata(hObject,handles);
    pause(0.2);
    handles = guidata(hObject);
end

% --- Executes on button press in Layer1.
function Layer1_Callback(hObject, eventdata, handles)
% hObject    handle to Layer1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Layer1
set(handles.Layer2,'Value',0);
set(handles.Layer3,'Value',0);
handles.I0 = handles.I(:,:,1);
imagesc(handles.I0);
if ~isempty(handles.Rect0)
    handles.Rect = rectangle('Position',handles.Rect0);
    set(handles.Rect,'EdgeColor','r','LineWidth',3);
end
if ~isempty(handles.tempRect1)
    handles.Rect2 = rectangle('Position',handles.tempRect1);
    set(handles.Rect2,'EdgeColor','r','LineWidth',3);
end
guidata(hObject,handles);

% --- Executes on button press in Layer2.
function Layer2_Callback(hObject, eventdata, handles)
% hObject    handle to Layer2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Layer2
set(handles.Layer1,'Value',0);
set(handles.Layer3,'Value',0);
handles.I0 = handles.I(:,:,2);
imagesc(handles.I0);
if ~isempty(handles.Rect0)
    handles.Rect = rectangle('Position',handles.Rect0);
    set(handles.Rect,'EdgeColor','r','LineWidth',3);
end
if ~isempty(handles.tempRect1)
    handles.Rect2 = rectangle('Position',handles.tempRect1);
    set(handles.Rect2,'EdgeColor','r','LineWidth',3);
end
guidata(hObject,handles);

% --- Executes on button press in Layer3.
function Layer3_Callback(hObject, eventdata, handles)
% hObject    handle to Layer3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Layer3
set(handles.Layer2,'Value',0);
set(handles.Layer1,'Value',0);
if size(handles.I,3)>2
    handles.I0 = handles.I(:,:,3);
    imagesc(handles.I0);
    if ~isempty(handles.Rect0)
        handles.Rect = rectangle('Position',handles.Rect0);
        set(handles.Rect,'EdgeColor','r','LineWidth',3);
    end
    if ~isempty(handles.tempRect1)
        handles.Rect2 = rectangle('Position',handles.tempRect1);
        set(handles.Rect2,'EdgeColor','r','LineWidth',3);
    end
else
    disp('image does not have a 3''rd layer');
    imagesc([nan nan nan; nan nan nan]);
end
guidata(hObject,handles);

% --- Executes on button press in LockRect2.
function LockRect2_Callback(hObject, eventdata, handles)
% hObject    handle to LockRect2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of LockRect2
handles = guidata(hObject);
if get(hObject,'Value')
    handles.tempRect1 = getPosition(handles.tempRect);
    handles.Rect2 = rectangle('Position',handles.tempRect1);
    set(handles.Rect2,'EdgeColor','r','LineWidth',3)
    guidata(hObject,handles);
end


% --- Executes on button press in countCells.
function countCells_Callback(hObject, eventdata, handles)
% hObject    handle to countCells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of countCells
% restrict image to boxes

if handles.countAll.Value % count cells for the entire image
    temp = adapthisteq(handles.I0);
    temp = double(temp);
    temp(temp==0)=nan; % 2019-11-05 AndyP
    thrVal0 = nanmean(temp(:))+handles.thrVal1*nanstd(temp(:));
    Ithr = temp > thrVal0;
    Ithr = medfilt2(Ithr,[3,3]);
    props = {'Area', 'PixelIdxList', 'Centroid', 'Eccentricity', 'Perimeter', 'Solidity'};
    BW = regionprops(Ithr,props);
    Cells = nan(length(BW),1);
    for iD=1:length(BW)
        if BW(iD).Area<handles.minSizeVal
            Cells(iD)=0;
        else
            Cells(iD)=1;
        end
    end
    Area = [];
    Centroid = [];
    Eccentricity = [];
    Solidity = [];
    for iD=1:length(BW)
        if Cells(iD)
            Area = cat(1,Area,BW(iD).Area);
            Centroid = cat(1,Centroid,BW(iD).Centroid);
            Eccentricity = cat(1,Eccentricity, BW(iD).Eccentricity);
            Solidity = cat(1,Solidity, BW(iD).Solidity);
        end
    end
    
    handles.areas1 = Area;
    handles.cell_pos1 = Centroid;
    handles.eccentricity1 = Eccentricity;
    handles.solidity1 = Solidity;
    
    figure(1);
    imagesc(handles.I0);
    hold on;
    plot(handles.cell_pos1(:,1),handles.cell_pos1(:,2),'r.','markersize',5);
    
else % count cells in rectangles
    pL = handles.Rect.Position;
    pR = handles.Rect2.Position;
    I1 = handles.I0(pL(2):pL(2)+pL(4),pL(1):pL(1)+pL(3));
    I2 = handles.I0(pR(2):pR(2)+pR(4),pR(1):pR(1)+pR(3));
    
    temp = double(I1);
    thrVal0 = nanmean(temp(:))+handles.thrVal1*nanstd(temp(:));
    Ithr = temp > thrVal0;
    props = {'Area', 'PixelIdxList', 'Centroid', 'Eccentricity', 'Perimeter', 'Solidity'};
    BW = regionprops(Ithr,props);
    Cells = nan(length(BW),1);
    for iD=1:length(BW)
        if BW(iD).Area<handles.minSizeVal
            Cells(iD)=0;
        else
            Cells(iD)=1;
        end
    end
    Area = [];
    Centroid = [];
    Eccentricity = [];
    Solidity = [];
    for iD=1:length(BW)
        if Cells(iD)
            Area = cat(1,Area,BW(iD).Area);
            Centroid = cat(1,Centroid,BW(iD).Centroid);
            Eccentricity = cat(1,Eccentricity, BW(iD).Eccentricity);
            Solidity = cat(1,Solidity, BW(iD).Solidity);
        end
    end
    
    handles.areas1 = Area;
    handles.cell_pos1 = Centroid;
    handles.eccentricity1 = Eccentricity;
    handles.solidity1 = Solidity;
    
    temp = double(I2);
    thrVal0 = nanmean(temp(:))+handles.thrVal1*nanstd(temp(:));
    Ithr = temp > thrVal0;
    props = {'Area', 'PixelIdxList', 'Centroid', 'Eccentricity', 'Perimeter', 'Solidity'};
    BW = regionprops(Ithr,props);
    Cells = nan(length(BW),1);
    for iD=1:length(BW)
        if BW(iD).Area<handles.minSizeVal
            Cells(iD)=0;
        else
            Cells(iD)=1;
        end
    end
    Area = [];
    Centroid = [];
    Eccentricity = [];
    Solidity = [];
    for iD=1:length(BW)
        if Cells(iD)
            Area = cat(1,Area,BW(iD).Area);
            Centroid = cat(1,Centroid,BW(iD).Centroid);
            Eccentricity = cat(1,Eccentricity, BW(iD).Eccentricity);
            Solidity = cat(1,Solidity, BW(iD).Solidity);
        end
    end
    
    handles.areas2 = Area;
    handles.cell_pos2 = Centroid;
    handles.eccentricity2 = Eccentricity;
    handles.solidity2 = Solidity;
    
    figure(1);
    imagesc(I1);
    hold on;
    plot(handles.cell_pos1(:,1),handles.cell_pos1(:,2),'r.','markersize',5);
    
    figure(2);
    imagesc(I2);
    hold on;
    plot(handles.cell_pos2(:,1),handles.cell_pos2(:,2),'r.','markersize',5);
end

%[handles.areas1,handles.cell_pos1]=countCells_2018_AP(I1,handles.sizeLims,handles.thrVals);

% figure(1);
% imagesc(I1);
% hold on;
% plot(handles.cell_pos1{1}(:,1),handles.cell_pos1{1}(:,2),'r.','markersize',5);
% 
% [handles.areas2,handles.cell_pos2]=countCells_2018_AP(I2,handles.sizeLims,handles.thrVals);
% 
% figure(2);
% imagesc(I2);
% hold on;
% plot(handles.cell_pos2{1}(:,1),handles.cell_pos2{1}(:,2),'r.','markersize',5);
% 
handles.countCells.Value=0;
guidata(hObject,handles);



function Thr1_Callback(hObject, eventdata, handles)
% hObject    handle to Thr1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Thr1 as text
%        str2double(get(hObject,'String')) returns contents of Thr1 as a double
handles.thrVal1 = str2double(get(handles.Thr1,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function Thr1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Thr1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function minSize_Callback(hObject, eventdata, handles)
% hObject    handle to minSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minSize as text
%        str2double(get(hObject,'String')) returns contents of minSize as a double
handles.minSizeVal = str2double(get(handles.minSize,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function minSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in countAll.
function countAll_Callback(hObject, eventdata, handles)
% hObject    handle to countAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of countAll
set(handles.CopyRect,'visible','off');
set(handles.drawRect,'visible','off');
set(handles.LockRect2,'visible','off');


% --- Executes on button press in memStain.
function memStain_Callback(hObject, eventdata, handles)
% hObject    handle to memStain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%handles.I0 = adapthisteq(handles.I0);
temp = handles.I0;
temp = adapthisteq(temp);
temp = double(temp);
thrVal0 = nanmean(temp(:))+handles.thrVal1*nanstd(temp(:));
temp = temp > thrVal0;
temp = medfilt2(temp,[3,3]);
Ithr = watershed(temp);
[handles.cell_pos1,handles.Radius] = imfindcircles(temp,[7 10],'Sensitivity',handles.SensVal,'EdgeThreshold',handles.EdgeThrVal,'Method','TwoStage');
figure;
imagesc(handles.I0);
hold on;
plot(handles.cell_pos1(:,1),handles.cell_pos1(:,2),'r.');
viscircles(handles.cell_pos1,handles.Radius);
% Hint: get(hObject,'Value') returns toggle state of memStain
% props = {'Area', 'PixelIdxList', 'Centroid', 'Eccentricity', 'Perimeter', 'Solidity'};
% BW = regionprops(Ithr,props);
% Cells = nan(length(BW),1);
% for iD=1:length(BW)
%     if BW(iD).Area<handles.minSizeVal
%         Cells(iD)=0;
%     else
%         Cells(iD)=1;
%     end
% end
% Area = [];
% Centroid = [];
% Eccentricity = [];
% Solidity = [];
% for iD=1:length(BW)
%     if Cells(iD)
%         Area = cat(1,Area,BW(iD).Area);
%         Centroid = cat(1,Centroid,BW(iD).Centroid);
%         Eccentricity = cat(1,Eccentricity, BW(iD).Eccentricity);
%         Solidity = cat(1,Solidity, BW(iD).Solidity);
%     end
% end
% 
% handles.areas1 = Area;
% handles.cell_pos1 = Centroid;
% handles.eccentricity1 = Eccentricity;
% handles.solidity1 = Solidity;
% 
% figure(1);
% imagesc(handles.I0);
% hold on;
% plot(handles.cell_pos1(:,1),handles.cell_pos1(:,2),'r.','markersize',5);

guidata(hObject,handles);

function Sensitivity_Callback(hObject, eventdata, handles)
% hObject    handle to Sensitivity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Sensitivity as text
%        str2double(get(hObject,'String')) returns contents of Sensitivity as a double
handles.SensVal = str2double(get(handles.Sensitivity,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function Sensitivity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Sensitivity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EdgeThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to EdgeThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EdgeThreshold as text
%        str2double(get(hObject,'String')) returns contents of EdgeThreshold as a double
handles.EdgeThrVal = str2double(get(handles.EdgeThreshold,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function EdgeThreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EdgeThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



