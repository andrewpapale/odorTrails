function varargout = getTrail_GUI3(varargin)
% 2017-06-13 AndyP
% getTrail_GUI3;
% used to extract odor trail from median image of optimouse output
% INPUT
% requires one input 1) arena_data.MedianImage, the median image from a
% video showing an odor trail
% OUTPUT
% saves *-Odor-Trail.mat with the coordinates (xT,yT) of the odor trail and
% a binary image of the odor trail outline (data) and filled region
% (filled).
% METHOD
% 1) Uses Canny edge detection and median filtering to extract trail based on user-set parameters from the GUI, 
% 2) uses a filling algorithm that connects the edges of the outline
% together and then fills in the shape of the odor trail.
% GETTRAIL_GUI3 MATLAB code for getTrail_GUI3.fig
%      GETTRAIL_GUI3, by itself, creates a new GETTRAIL_GUI3 or raises the existing
%      singleton*.
%
%      H = GETTRAIL_GUI3 returns the handle to a new GETTRAIL_GUI3 or the handle to
%      the existing singleton*.
%
%      GETTRAIL_GUI3('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GETTRAIL_GUI3.M with the given input arguments.
%
%      GETTRAIL_GUI3('Property','Value',...) creates a new GETTRAIL_GUI3 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before getTrail_GUI3_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to getTrail_GUI3_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help getTrail_GUI3

% Last Modified by GUIDE v2.5 14-Jun-2017 14:05:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @getTrail_GUI3_OpeningFcn, ...
    'gui_OutputFcn',  @getTrail_GUI3_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
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


% --- Executes just before getTrail_GUI3 is made visible.
function getTrail_GUI3_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to getTrail_GUI3 (see VARARGIN)

% Choose default command line output for getTrail_GUI3
handles.output = hObject;

[handles.FileName,~,~] = uigetfile('.jpg');
mI = double(imread(handles.FileName));
white = (mI(:,:,1)==255 & mI(:,:,2)==255 & mI(:,:,3)==255);
mI = squeeze((mI(:,:,1)-mI(:,:,2)));
mI(white)=255;
%mI = double(arena_data.MedianImage);
handles.current_data = mI;
imagesc(handles.current_data);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes getTrail_GUI3 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = getTrail_GUI3_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
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

nbrh = getRadioButton(handles);

data = getData_getTrailGUI(handles,nbrh);

tempStr = strsplit(handles.FileName,'_');
saveStr = strcat(tempStr{1},'_',tempStr{2},'_',tempStr{3},'_',tempStr{4},'_',tempStr{5},'-Odor-Trail.mat');
[xT,yT]=find(data==1); %#ok<ASGLU>
fprintf('Saving %s \n',saveStr);
save(saveStr,'data','xT','yT','minThr','maxThr','sigma','nbrh','medianFilter','Thr','edgeVal','edgeDetection','Filling','ThresholdData');
close all force



% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.radiobutton1.Value
    handles.radiobutton7.Value = 0;
    handles.radiobutton5.Value = 0;
    handles.radiobutton4.Value = 0;
    handles.radiobutton3.Value = 0;
    handles.radiobutton2.Value = 0;
    handles.radiobutton6.Value = 0;
end
nbrh = getRadioButton(handles);
data = getData_getTrailGUI(handles,nbrh);
imagesc(data);
% Hint: get(hObject,'Value') returns toggle state of radiobutton1


% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.radiobutton2.Value
    handles.radiobutton7.Value = 0;
    handles.radiobutton5.Value = 0;
    handles.radiobutton4.Value = 0;
    handles.radiobutton3.Value = 0;
    handles.radiobutton6.Value = 0;
    handles.radiobutton1.Value = 0;
end

nbrh = getRadioButton(handles);
data = getData_getTrailGUI(handles,nbrh);
imagesc(data);

% Hint: get(hObject,'Value') returns toggle state of radiobutton2


% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.radiobutton3.Value
    handles.radiobutton7.Value = 0;
    handles.radiobutton5.Value = 0;
    handles.radiobutton4.Value = 0;
    handles.radiobutton6.Value = 0;
    handles.radiobutton2.Value = 0;
    handles.radiobutton1.Value = 0;
end


nbrh = getRadioButton(handles);
data = getData_getTrailGUI(handles,nbrh);
imagesc(data);
% Hint: get(hObject,'Value') returns toggle state of radiobutton3


% --- Executes on button press in radiobutton4.
function radiobutton4_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.radiobutton4.Value
    handles.radiobutton7.Value = 0;
    handles.radiobutton5.Value = 0;
    handles.radiobutton6.Value = 0;
    handles.radiobutton3.Value = 0;
    handles.radiobutton2.Value = 0;
    handles.radiobutton1.Value = 0;
end

nbrh = getRadioButton(handles);
data = getData_getTrailGUI(handles,nbrh);
imagesc(data);

% Hint: get(hObject,'Value') returns toggle state of radiobutton4


% --- Executes on button press in radiobutton5.
function radiobutton5_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.radiobutton5.Value
    handles.radiobutton7.Value = 0;
    handles.radiobutton6.Value = 0;
    handles.radiobutton4.Value = 0;
    handles.radiobutton3.Value = 0;
    handles.radiobutton2.Value = 0;
    handles.radiobutton1.Value = 0;
end

nbrh = getRadioButton(handles);
data = getData_getTrailGUI(handles,nbrh);
imagesc(data);

% Hint: get(hObject,'Value') returns toggle state of radiobutton5


% --- Executes on button press in radiobutton6.
function radiobutton6_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.radiobutton6.Value
    handles.radiobutton7.Value = 0;
    handles.radiobutton5.Value = 0;
    handles.radiobutton4.Value = 0;
    handles.radiobutton3.Value = 0;
    handles.radiobutton2.Value = 0;
    handles.radiobutton1.Value = 0;
end

nbrh = getRadioButton(handles);
data = getData_getTrailGUI(handles,nbrh);
imagesc(data);
% Hint: get(hObject,'Value') returns toggle state of radiobutton6


% --- Executes on button press in radiobutton7.
function radiobutton7_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.radiobutton7.Value
    handles.radiobutton6.Value = 0;
    handles.radiobutton5.Value = 0;
    handles.radiobutton4.Value = 0;
    handles.radiobutton3.Value = 0;
    handles.radiobutton2.Value = 0;
    handles.radiobutton1.Value = 0;
end

nbrh = getRadioButton(handles);
data = getData_getTrailGUI(handles,nbrh);
imagesc(data);
% Hint: get(hObject,'Value') returns toggle state of radiobutton7


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.checkbox1,'Value');
    set(handles.checkbox4,'Value',0);
    set(handles.checkbox4,'Enable','off');
else
    set(handles.checkbox4,'Enable','on');
end
if ~get(handles.checkbox1,'Value') && ~get(handles.checkbox4,'Value');
    set(handles.checkbox3,'Enable','off');
else
    set(handles.checkbox3,'Enable','on');
end
nbrh = getRadioButton(handles);
data = getData_getTrailGUI(handles,nbrh);
imagesc(data);
% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

nbrh = getRadioButton(handles);
data = getData_getTrailGUI(handles,nbrh);
imagesc(data);
% Hint: get(hObject,'Value') returns toggle state of checkbox2


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
set(handles.slider1,'max',get(handles.slider3,'Value')-0.01);
set(handles.edit1,'String',get(handles.slider1,'Value'));
set(handles.edit2,'String',get(handles.slider3,'Value'));
set(handles.edit3,'String',get(handles.slider4,'Value'));


nbrh = getRadioButton(handles);
data = getData_getTrailGUI(handles,nbrh);
imagesc(data);

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.slider1,'max',get(handles.slider3,'Value')-0.01);
set(handles.edit1,'String',get(handles.slider1,'Value'));
set(handles.edit2,'String',get(handles.slider3,'Value'));
set(handles.edit3,'String',get(handles.slider4,'Value'));

nbrh = getRadioButton(handles);
data = getData_getTrailGUI(handles,nbrh);
imagesc(data);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider4_Callback(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.slider1,'max',get(handles.slider3,'Value')-0.01);
set(handles.edit1,'String',get(handles.slider1,'Value'));
set(handles.edit2,'String',get(handles.slider3,'Value'));
set(handles.edit3,'String',get(handles.slider4,'Value'));


nbrh = getRadioButton(handles);
data = getData_getTrailGUI(handles,nbrh);
imagesc(data);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
nbrh = getRadioButton(handles);
data = getData_getTrailGUI(handles,nbrh);
imagesc(data);
% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

nbrh = getRadioButton(handles);
data = getData_getTrailGUI(handles,nbrh);
imagesc(data);
% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

nbrh = getRadioButton(handles);
data = getData_getTrailGUI(handles,nbrh);
imagesc(data);
% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double

nbrh = getRadioButton(handles);
data = getData_getTrailGUI(handles);
imagesc(data);
% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
nbrh = getRadioButton(handles);
data = getData_getTrailGUI(handles,nbrh);
imagesc(data);


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(handles.checkbox4,'Value');
    set(handles.checkbox1,'Value',0);
    set(handles.checkbox1,'Enable','off');
else
    set(handles.checkbox1,'Enable','on');
end
if ~get(handles.checkbox1,'Value') && ~get(handles.checkbox4,'Value');
    set(handles.checkbox3,'Enable','off');
else
    set(handles.checkbox3,'Enable','on');
end
nbrh = getRadioButton(handles);
data = getData_getTrailGUI(handles,nbrh);
imagesc(data);


% Hint: get(hObject,'Value') returns toggle state of checkbox4


% --- Executes on slider movement.
function slider5_Callback(hObject, eventdata, handles)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
set(handles.edit5,'String',get(handles.slider5,'Value'));
nbrh = getRadioButton(handles);
data = getData_getTrailGUI(handles,nbrh);
imagesc(data);

% --- Executes during object creation, after setting all properties.
function slider5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double
nbrh = getRadioButton(handles);
data = getData_getTrailGUI(handles,nbrh);
imagesc(data);



% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider6_Callback(hObject, eventdata, handles)
% hObject    handle to slider6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
set(handles.edit6,'String',get(handles.slider6,'Value'));
nbrh = getRadioButton(handles);
data = getData_getTrailGUI(handles,nbrh);
imagesc(data);

% --- Executes during object creation, after setting all properties.
function slider6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.edit6,'String',get(handles.slider6,'Value'));
nbrh = getRadioButton(handles);
data = getData_getTrailGUI(handles,nbrh);
imagesc(data);

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
