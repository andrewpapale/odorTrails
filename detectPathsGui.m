function varargout = detectPathsGui(varargin)
% DETECTPATHSGUI MATLAB code for detectPathsGui.fig
%      DETECTPATHSGUI, by itself, creates a new DETECTPATHSGUI or raises the existing
%      singleton*.
%
%      H = DETECTPATHSGUI returns the handle to a new DETECTPATHSGUI or the handle to
%      the existing singleton*.
%
%      DETECTPATHSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DETECTPATHSGUI.M with the given input arguments.
%
%      DETECTPATHSGUI('Property','Value',...) creates a new DETECTPATHSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before detectPathsGui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to detectPathsGui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help detectPathsGui

% Last Modified by GUIDE v2.5 18-Mar-2014 16:22:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @detectPathsGui_OpeningFcn, ...
                   'gui_OutputFcn',  @detectPathsGui_OutputFcn, ...
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


% --- Executes just before detectPathsGui is made visible.
function detectPathsGui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to detectPathsGui (see VARARGIN)

% Choose default command line output for detectPathsGui
handles.output = hObject;
global VIDEO_ROOT;
handles.basePath = VIDEO_ROOT;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes detectPathsGui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = detectPathsGui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in load_button.
function load_button_Callback(hObject, eventdata, handles)
% hObject    handle to load_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%base_path = VIDEO_ROOT;
PathName = uigetdir(handles.basePath, 'Select the folder to process');
folders = textscan(PathName, '%s', 'Delimiter', filesep);
handles.folders = folders{1};
handles.exp_name = folders{end};
if iscell(handles.exp_name)
    handles.exp_name = handles.exp_name{end};
end
handles.vids = processVideoFolder(handles.exp_name, @MouseTrackerKF, 0);
handles.vidi = 1;

handles = setVideoIndex(hObject, handles, 1);
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in copy_button.
function copy_button_Callback(hObject, eventdata, handles)
% hObject    handle to copy_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.vidi > 1
    handles.vids(handles.vidi).paths = handles.vids(handles.vidi-1).paths;
    handles.vids(handles.vidi).save;
end
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in detect_button.
function detect_button_Callback(hObject, eventdata, handles)
% hObject    handle to detect_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.vids(handles.vidi).detectRefinePaths(1,1,1, handles.im_ax);
handles.vids(handles.vidi).save;
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in next_button.
function next_button_Callback(hObject, eventdata, handles)
% hObject    handle to next_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.vidi < length(handles.vids)
    handles = setVideoIndex(hObject, handles, handles.vidi + 1);
end
    % Update handles structure
guidata(hObject, handles);


function videoNum_text_Callback(hObject, eventdata, handles)
% hObject    handle to videoNum_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of videoNum_text as text
%        str2double(get(hObject,'String')) returns contents of videoNum_text as a double

handles = setVideoIndex(hObject, handles, str2double(get(hObject,'String')));
    % Update handles structure
guidata(hObject, handles);



function handles = setVideoIndex(hObject, handles, ind)
%function setVideoIndex(hObject, handles, ind)
%
% Called by several functions to set the video index and do proper updating of GUI

handles.vidi = ind;
set(handles.videoNum_text, 'String', num2str(ind));
axes(handles.im_ax);
pathIm = handles.vids(handles.vidi).plotPathsOnBg;
imshow(pathIm);
% Update handles structure
%guidata(hObject, handles);