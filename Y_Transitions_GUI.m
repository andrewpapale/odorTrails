function varargout = Y_Transitions_GUI(varargin)
% Y_TRANSITIONS_GUI MATLAB code for Y_Transitions_GUI.fig
%      Y_TRANSITIONS_GUI, by itself, creates a new Y_TRANSITIONS_GUI or raises the existing
%      singleton*.
%
%      H = Y_TRANSITIONS_GUI returns the handle to a new Y_TRANSITIONS_GUI or the handle to
%      the existing singleton*.
%
%      Y_TRANSITIONS_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in Y_TRANSITIONS_GUI.M with the given input arguments.
%
%      Y_TRANSITIONS_GUI('Property','Value',...) creates a new Y_TRANSITIONS_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Y_Transitions_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Y_Transitions_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Y_Transitions_GUI

% Last Modified by GUIDE v2.5 03-Oct-2018 13:21:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Y_Transitions_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @Y_Transitions_GUI_OutputFcn, ...
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


% --- Executes just before Y_Transitions_GUI is made visible.
function Y_Transitions_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Y_Transitions_GUI (see VARARGIN)

% Choose default command line output for Y_Transitions_GUI
handles.output = hObject;

[handles.FileName,~,~] = uigetfile;
load(handles.FileName);
tempStr = strsplit(handles.FileName,'_');
tempStr1 = strcat(tempStr{1},'_',tempStr{2},'_',tempStr{3},'_',tempStr{4},'_',tempStr{5},'-ConnectingPoint.mat');
load(tempStr1,'xC','yC');
handles.xC = xC;
handles.yC = yC;
tempStr2 = strcat(tempStr{1},'_',tempStr{2}(1:17),'-StartFrame.mat');
load(tempStr2);
tempStr3 = strcat(tempStr{1},'_',tempStr{2}(1:22),'_',tempStr{3},'_',tempStr{4},'_',tempStr{5},'-Odor-Trail.mat');
load(tempStr3,'xT','yT');
handles.yT = xT;
handles.xT = yT;
% calculate passes through center point <xC,yC>
[handles.x0,handles.y0,handles.nx0,handles.ny0,~,~,~,~,~] = Process_VT(position_results,startFrame,1);
dnTc = sqrt((handles.x0-handles.xC).^2+(handles.y0-handles.yC).^2)./11.2;
handles.incenter = find(dnTc < 10);
handles.passes = diff(handles.incenter)>1;
handles.nP = sum(handles.passes);
handles.iP = 1;

set(handles.PassCount,'String',sprintf('Pass %d/%d',handles.iP,handles.nP));
plot(handles.xT,handles.yT,'r.'); 
hold on;
plot(handles.x0,handles.y0,'k.');
k1 = 1;
k2 = find(handles.passes==1);
k2 = k2(handles.iP);
handles.P1 = scatter(handles.x0(handles.incenter(k1:k2)),handles.y0(handles.incenter(k1:k2)),20,1:sum(k2-k1)+1,'filled');
set(gca,'XLim',[handles.xC-200,handles.xC+200],'YLim',[handles.yC-200,handles.yC+200]);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Y_Transitions_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Y_Transitions_GUI_OutputFcn(hObject, eventdata, handles) 
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
close all force;

% --- Executes on button press in SaveData.
function SaveData_Callback(hObject, eventdata, handles)
% hObject    handle to SaveData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in NextTran.
function NextTran_Callback(hObject, eventdata, handles)
% hObject    handle to NextTran (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.iP~=handles.nP
    delete(handles.P1);
    handles.iP = handles.iP+1;
end
set(handles.PassCount,'String',sprintf('Pass %d/%d',handles.iP,handles.nP));
plot(handles.xT,handles.yT,'r.'); 
hold on;
plot(handles.x0,handles.y0,'k.');
k1 = find(handles.passes==1);
k2 = k1(handles.iP);
k1 = k1(handles.iP-1);
handles.P1 = scatter(handles.x0(handles.incenter(k1:k2)),handles.y0(handles.incenter(k1:k2)),20,1:sum(k2-k1)+1,'filled');
set(gca,'XLim',[handles.xC-200,handles.xC+200],'YLim',[handles.yC-200,handles.yC+200]);
guidata(hObject,handles);

% --- Executes on button press in PrevTran.
function PrevTran_Callback(hObject, eventdata, handles)
% hObject    handle to PrevTran (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.iP>1
    delete(handles.P1);
    handles.iP = handles.iP-1;
    set(handles.PassCount,'String',sprintf('Pass %d/%d',handles.iP,handles.nP));
end
plot(handles.xT,handles.yT,'r.'); 
hold on;
plot(handles.x0,handles.y0,'k.');
k1 = find(handles.passes==1);
k2 = k1(handles.iP);
k1 = k1(handles.iP-1);
handles.P1 = scatter(handles.x0(handles.incenter(k1:k2)),handles.y0(handles.incenter(k1:k2)),20,1:sum(k2-k1)+1,'filled');
set(gca,'XLim',[handles.xC-200,handles.xC+200],'YLim',[handles.yC-200,handles.yC+200]);
guidata(hObject,handles);

% --- Executes on button press in One2Two.
function One2Two_Callback(hObject, eventdata, handles)
% hObject    handle to One2Two (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of One2Two


% --- Executes on button press in Two2One.
function Two2One_Callback(hObject, eventdata, handles)
% hObject    handle to Two2One (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Two2One


% --- Executes on button press in One2Three.
function One2Three_Callback(hObject, eventdata, handles)
% hObject    handle to One2Three (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of One2Three


% --- Executes on button press in Three2One.
function Three2One_Callback(hObject, eventdata, handles)
% hObject    handle to Three2One (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Three2One


% --- Executes on button press in Two2Three.
function Two2Three_Callback(hObject, eventdata, handles)
% hObject    handle to Two2Three (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Two2Three


% --- Executes on button press in Three2Two.
function Three2Two_Callback(hObject, eventdata, handles)
% hObject    handle to Three2Two (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Three2Two


% --- Executes on button press in NoneofTheAbove.
function NoneofTheAbove_Callback(hObject, eventdata, handles)
% hObject    handle to NoneofTheAbove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of NoneofTheAbove



function PassCount_Callback(hObject, eventdata, handles)
% hObject    handle to PassCount (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PassCount as text
%        str2double(get(hObject,'String')) returns contents of PassCount as a double


% --- Executes during object creation, after setting all properties.
function PassCount_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PassCount (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
