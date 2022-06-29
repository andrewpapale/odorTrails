function varargout = Tenzin_GUI(varargin)
%Tenzin_GUI MATLAB code file for Tenzin_GUI.fig
%      Tenzin_GUI, by itself, creates a new TARA_GUI or raises the existing
%      singleton*.
%
%      H = Tenzin_GUI returns the handle to a new TARA_GUI or the handle to
%      the existing singleton*.
%
%      Tenzin_GUI('Property','Value',...) creates a new TARA_GUI using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to Tenzin_GUI_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      Tenzin_GUI('CALLBACK') and Tenzin_GUI('CALLBACK',hObject,...) call the
%      local function named CALLBACK in Tenzin_GUI.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Tenzin_GUI

% Last Modified by GUIDE v2.5 16-Jul-2019 12:35:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Tenzin_GUI_OpeningFcn, ...
    'gui_OutputFcn',  @Tenzin_GUI_OutputFcn, ...
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


% --- Executes just before Tenzin_GUI is made visible.
function Tenzin_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for Tenzin_GUI
handles.output = hObject;
[handles.FileName,~,~] = uigetfile;
handles.t = Tiff(handles.FileName,'r');
% handles.thrVals = [0.05 0.1 0.15];
% handles.sizeLims = [50 140];
handles.thrVal1 = 2.5;
handles.minSizeVal = 40;
handles.pixperum = 35/50;
handles.cellsCounted = 0;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Tenzin_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Tenzin_GUI_OutputFcn(hObject, eventdata, handles)
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


% --- Executes on button press in SaveData.
function SaveData_Callback(hObject, eventdata, handles)
% hObject    handle to SaveData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
areas1 = handles.areas1;
cell_pos1 = handles.cell_pos1;
eccentricity1 = handles.eccentricity1;
solidity1 = handles.solidity1;
threshold = handles.thrVal1;
minSize = handles.minSizeVal;

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
save(saveStr,'areas1','cell_pos1','areas2','cell_pos2','eccentricity1','eccentricity2','solidity1','solidity2','threshold','minSize');


% --- Executes on button press in countCells.
function countCells_Callback(hObject, eventdata, handles)
% hObject    handle to countCells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of countCells
% restrict image to boxes
nZ = handles.t.numberOfStrips;
Area = [];
Centroid = [];
Eccentricity = [];
Solidity = [];
Zmatrix = [];
for iZ=1:nZ
    iC = 1;
    handles.zP = iZ;
    handles.I0 = imread(handles.FileName,handles.zP);
    temp = double(handles.I0-nanmedian(handles.I0));
    thrVal0 = nanmean(temp(:))+handles.thrVal1*nanstd(temp(:));
    Ithr = temp > thrVal0;
    props = {'Area', 'PixelIdxList', 'Centroid', 'Eccentricity','Perimeter', 'Solidity'};
    BW = regionprops(Ithr,props);
    Cells = nan(length(BW),1);
    for iD=1:length(BW)
        if BW(iD).Area<handles.minSizeVal
            Cells(iD)=0;
        else
            Cells(iD)=1;
        end
    end
    for iD=1:length(BW)
        if Cells(iD)
            Zmatrix{iZ,iC} = iZ;
            Area{iZ,iC} = BW(iD).Area;
            Centroid{iZ,iC} = BW(iD).Centroid;
            Eccentricity{iZ,iC} = BW(iD).Eccentricity;
            Solidity{iZ,iC} = BW(iD).Solidity;
            iC = iC+1;
        end
    end
    imagesc(handles.I0);
    hold on;
    if ~isempty(Centroid)
        nCz = sum(Cells);
        for iC0=1:nCz
            C0 = Centroid{iZ,iC0};
            plot(C0(:,1),C0(:,2),'rx','markersize',10);
            hold on;
        end
    end
    handles = guidata(hObject);
    guidata(hObject,handles);
    pause(0.1);
end
% go through and eliminate doubles (within
mZ = [];
a0 = [];
cx0 = [];
cy0 = [];
e0 = [];
s0 = [];
nZ = size(Zmatrix,1);
for iZ=1:nZ
    nC0 = size(Zmatrix,2);
    if ~isempty(nC0)
        for iC=1:nC0
            if ~isempty(Zmatrix{iZ,iC})
                mZ = cat(1,mZ,Zmatrix{iZ,iC});
                a0 = cat(1,a0,Area{iZ,iC});
                cx0 = cat(1,cx0,Centroid{iZ,iC}(1,1));
                cy0 = cat(1,cy0,Centroid{iZ,iC}(1,2));
                e0 = cat(1,e0,Eccentricity{iZ,iC});
                s0 = cat(1,s0,Solidity{iZ,iC});
            end
        end
    end
end
% construct distance between all points in each plane
remove = zeros(size(mZ));
iC = 1;
for iZ=1:nZ
    k = find(mZ==iZ);
    if ~isempty(k)
        nP = length(k);
        nQ = nP;
        Dxy = nan(nP,nQ);
        cx1 = cx0(k);
        cy1 = cy0(k);
        for iP=1:nP
            for iQ=1:nQ-iP
                if iQ~=iP
                    Dxy(iP,iQ)=sqrt(((cx1(iP)-cx1(iQ)).^2)./handles.pixperum+((cy1(iP)-cy1(iQ)).^2)./handles.pixperum);
                end
            end
            if min(Dxy(iP,:))<15
                remove(iC)=1;
                iC = iC+1;
            end
        end
    end
end
a0(remove==1)=[];
cx0(remove==1)=[];
cy0(remove==1)=[];
e0(remove==1)=[];
s0(remove==1)=[];
mZ(remove==1)=[];

% already threw away distance < threshold in plane z
% for each point z1 in plane z, determine distances < threshold in plane
% z+1
% consider closest point z2 in z+1 < threshold 'same cell' as point in z
% discard points z3...zN in z+1 < threshold but > distance of z2 as noise/axon labeling/etc
% 'keep' point z2, and repeat analysis for plane z+2...z+N

cI = nan(size(mZ)); % cell identity
iC = 0;
minZ = find(min(mZ)); % plane z
for iZ=minZ:nZ-1
    nPiz = find(mZ==iZ);
    nPiz1= find(mZ==iZ+1);
    if ~isempty(nPiz) && ~isempty(nPiz1)
        % compute distances
        nP = length(nPiz);
        nQ = length(nPiz1);
        Dxy = nan(nP,nQ);
        for iP=1:nP
            for iQ=1:nQ-iP
                if iP~=iQ
                    Dxy(iP,iQ)=sqrt((nx0(nPiz(iP))-nx0(nPiz1(iQ))).^2./handles.pixperum+(ny0(nPiz(iP))-ny0(nPiz1(iQ))).^2./handles.pixperum);
                end
            end
        end
        sameFlag = zeros(nP,1);
        while any(min(Dxy(:))) < 30
            [minDxy,NN] = min(Dxy,[],2); %
            for iP=1:nP
                if minDxy(iP) < 30 && ~sameFlag(iP) % classify as same
                    if isnan(nPiz1(NN)) % if not classified as cell  yet
                        cI(nPiz1(NN))=iC;
                    end
                    if isnan(nPiz(iP))
                        cI(nPiz(iP))=iC;
                    end
                    Dxy(iP,NN) = nan;
                    sameFlag(iP)=1;
                elseif sameFlag(iP) % noise point, remove
                    Dxy(iP,NN)=nan;
                else % > threshold, different cell
                    if isnan(nPiz(iP))
                        iC = iC+1;
                        cI(nPiz(iP))=iC;
                    end
                    if isnan(nPiz1(NN))
                        iC = iC+1;
                        cI(nPiz1(NN))=iC;
                    end
                    
                end
            end
        end
    end
end

    


handles.zstack = mZ;
handles.areas = a0;
handles.cell_pos = [cx0,cy0];
handles.eccentricity = e0;
handles.solidity = s0;
handles.cellsCounted = 1;
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

% --- Executes on slider movement.
function zSlider_Callback(hObject, eventdata, handles)
% hObject    handle to zSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles = guidata(hObject);
handles.zP = round(get(hObject, 'Value'));
handles.I0 = imread(handles.FileName,handles.zP);
imagesc(handles.I0);
if handles.cellsCounted
    if any(handles.zstack==handles.zP)
        hold on;
        ix = find(handles.zstack==handles.zP);
        plot(handles.cell_pos(ix,1),handles.cell_pos(ix,2),'r.','markersize',20);
    end
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function zSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject, 'Max', 100, 'Min', 1,'SliderStep',[1/99 0.1],'Value',1);

%handles.zP = round(get(hObject, 'Value'));


function zStackNum_Callback(hObject, eventdata, handles)
% hObject    handle to zStackNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zStackNum as text
%        str2double(get(hObject,'String')) returns contents of zStackNum as a double
handles = guidata(hObject);
str = get(handles.zStackNum,'String');
try
    handles.zP = round(str2double(str));
catch
    warning('this value must be a string');
end
set(handles.zSlider,'Value',handles.zP);
handles.I0 = imread(handles.FileName,handles.zP);
imagesc(handles.I0);
if handles.cellsCounted
    if any(handles.zstack==handles.zP)
        hold on;
        ix = find(handles.zstack==handles.zP);
        plot(handles.cell_pos(ix,1),handles.cell_pos(ix,2),'r.','markersize',20);
    end
end
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function zStackNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zStackNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

