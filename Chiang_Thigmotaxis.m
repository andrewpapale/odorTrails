function varargout = Chiang_Thigmotaxis(varargin)
% CHIANG_THIGMOTAXIS MATLAB code for Chiang_Thigmotaxis.fig
% 2018-07-21 AndyP
% 2018-10-11 AndyP, added single (left) box functionality, fixed pointer to
% edge val (was fixed at 5cm)
% GUI for Michael Chiang's analysis of data
% General Procedure:
% 1) Run "Chiang_Thigmotaxis" in the matlab command line.  This will promot
% the user to open a positions file (optimouse output).
% 2) Click on "Draw Left Rectangle" button.  This will prompt the user to
% draw a rectangle on the left side of the chamber
% 3) Click on "Draw Div Line" button.  This will prompt the user to click on two
% points, defining a line dividing the center of the chamber.
% 4) Click on "Draw Right Rectangle" button.  This will pop up a rectangle
% with equal area.  Drag the rectangle to the right side of the chamber.
% 5) When the right rectangle is properly positioned, click the checkbox "Rect OK"
% 6) Enter in a value for the edge threshold (in cm).  This will determine
% at what distance from the edge is considered "near the edge" versus "in
% the center of the chamber."  I've been using "5" cm.
% 7) Click the checkbox "Compute Edges".  This will plot the center and
% edge trajectories in different colors.
% 8) Click the checkbox "Compute Crossings".  This will compute the times
% the mouse crosses from LtoR or RtoL of the chamber.
% 9) If you want to start over, click on "Reset".
% 10) If everything is OK, click, "Save".

% Saved variables
% 'EdgeThrVal': The value you enter in step 6), determines edge threshold.
% 'Area': Area of the rectangles.
% 'RectL','RectR','Line': Data on the position of the rectangles and line.
% 'edgeL','edgeR','centL','centR': logicals to restrict position data.
% Same size as <x,y> data.  edgeL is left edge, etc.
% 'LtoL','LtoR','RtoR','RtoL','inL','inR': logicals to restrict position
% data.  Same size as <x,y> data.  inL are the samples in the left side,
% etc.
% 'nLR','nRL': the number of crossings from RtoL (nRL) and LtoR (nLR)
% 'tLR','tRL','tR','tL','tRedge','tLedge','tRcent','tLcent':  times (in s)
% when the mouse crosses from LtoR (tLR) RtoL (tRL), the total time in R
% (tR) and the total time in L (tL).

% Edit the above text to modify the response to help Chiang_Thigmotaxis

% Last Modified by GUIDE v2.5 09-Oct-2018 10:50:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Chiang_Thigmotaxis_OpeningFcn, ...
    'gui_OutputFcn',  @Chiang_Thigmotaxis_OutputFcn, ...
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


% --- Executes just before Chiang_Thigmotaxis is made visible.
function Chiang_Thigmotaxis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Chiang_Thigmotaxis (see VARARGIN)

% Choose default command line output for Chiang_Thigmotaxis
handles.output = hObject;
[handles.FileName,~,~] = uigetfile;
load(handles.FileName);
handles.mI = double(arena_data.MedianImage);
handles.x = position_results.mouseCOM(:,1);
handles.y = position_results.mouseCOM(:,2);
handles.pixpermm = arena_data.pixels_per_mm;
handles.RectL = [];
handles.RectR = [];
handles.Line = [];
handles.tempRect = [];
handles.LtoL = [];
handles.LtoR = [];
handles.RtoR = [];
handles.RtoL = [];
handles.inL = [];
handles.inR = [];
handles.nLR = [];
handles.nRL = [];
handles.tLR = [];
handles.tRL = [];
handles.tL_binned = [];
handles.tR_binned = [];
handles.edgeL_binned = [];
handles.edgeR_binned = [];
handles.centL_binned = [];
handles.centR_binned = [];
handles.ExclusionZone = [];
imagesc(handles.mI);
colormap gray;
hold on;
plot(handles.x,handles.y,'k.-');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Chiang_Thigmotaxis wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Chiang_Thigmotaxis_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Save.
function Save_Callback(hObject, eventdata, handles)
% hObject    handle to Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


tempStr = strsplit(handles.FileName,'_');
saveStr = strcat(tempStr{1},'_','positions','_','AnalyzedData.mat');

EdgeThrVal = handles.EdgeThrVal;
Area = handles.Area;
RectL = handles.RectL;
RectR = handles.RectR;
Line = handles.Line;
edgeL = handles.edgeL;
edgeR = handles.edgeR;
centL = handles.centL;
centR = handles.centR;
LtoL = handles.LtoL;
LtoR = handles.LtoR;
RtoR = handles.RtoR;
RtoL = handles.RtoL;
inL = handles.inL;
inR = handles.inR;
nLR = handles.nLR;
nRL = handles.nRL;
tLR = handles.tLR;
tRL = handles.tRL;
if ~isempty(RtoL)
    tR = sum(inR & ~RtoL & ~LtoR)/30;
    tL = sum(inL & ~RtoL & ~LtoR)/30;
else
    tR = 0;
    tL = sum(inL)/30;
end
if ~isempty(inR)
    if ~isempty(RtoL)
        tRedge = sum(edgeR & inR & ~RtoL & ~LtoR)/30;
        tRcent = sum(centR & inR & ~RtoL & ~LtoR)/30;
    else
        tRedge = sum(edgeR & inR)/30;
        tRcent = sum(centR & inR)/30;
    end
else
    tRedge = [];
    tRcent = [];
end
if ~isempty(inL)
    if ~isempty(RtoL)
        tLedge = sum(edgeL & inL & ~RtoL & ~LtoR)/30;
        tLcent = sum(centL & inL & ~RtoL & ~LtoR)/30;
    else
        tLedge = sum(edgeL & inL)/30;
        tLcent = sum(centL & inL)/30;
    end
else
    tLedge = [];
    tLcent = [];
end
tL_binned = handles.tL_binned;
tR_binned = handles.tR_binned;
edgeL_binned = handles.edgeL_binned;
edgeR_binned = handles.edgeR_binned;
centL_binned = handles.centL_binned;
centR_binned = handles.centR_binned;

save(saveStr,'EdgeThrVal','Area','RectL','RectR','Line','edgeL','edgeR','centL','centR',...
    'LtoL','LtoR','RtoR','RtoL','inL','inR','nLR','nRL','tLR','tRL','tR','tL','tRedge',...
    'tLedge','tRcent','tLcent','tL_binned','tR_binned','edgeL_binned','edgeR_binned',...
    'centL_binned','centR_binned');


% --- Executes on button press in Exit.
function Exit_Callback(hObject, eventdata, handles)
% hObject    handle to Exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.checkbox1.Value = 0;
close all force

% --- Executes on button press in RectL.
function RectL_Callback(hObject, eventdata, handles)
% hObject    handle to RectL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.RectR)
    handles.RectL = getrect;
    handles.RectL = rectangle('Position',handles.RectL);
    set(handles.RectL,'EdgeColor','r','LineWidth',3);
    guidata(hObject,handles);
else
    handles.tempRect = imrect(gca,handles.RectR.Position);
    handles.LorR = 'L';
    while ~handles.checkbox1.Value
        guidata(hObject,handles);
        pause(0.2);
        handles = guidata(hObject);
    end
    guidata(hObject,handles);
end

handles.Area = handles.RectL.Position(3)*handles.RectL.Position(4)/(handles.pixpermm*10).^2;
set(handles.AreaVal,'string',mat2str(round(handles.Area)));

guidata(hObject,handles);

% --- Executes on button press in RectR.
function RectR_Callback(hObject, eventdata, handles)
% hObject    handle to RectR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.RectL)
    handles.RectR = getrect;
    handles.RectR = rectangle('Position',handles.RectR);
    set(handles.RectR,'EdgeColor','r','LineWidth',3);
    guidata(hObject,handles);
else
    handles.tempRect = imrect(gca,handles.RectL.Position);
    handles.LorR = 'R';
    while ~handles.checkbox1.Value
        guidata(hObject,handles);
        pause(0.2);
        handles = guidata(hObject);
    end
end

handles.Area = handles.RectR.Position(3)*handles.RectR.Position(4)/(handles.pixpermm*10).^2;
set(handles.AreaVal,'string',mat2str(round(handles.Area)));

guidata(hObject,handles);

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = ginput(2);
handles.Line = line(data(:,1),data(:,2),'color','r','linewidth',3);

guidata(hObject,handles);

% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
handles = guidata(hObject);
if get(hObject,'Value')
    if ~isempty(handles.tempRect)
        handles.tempRect1 = getPosition(handles.tempRect);
        if strcmp(handles.LorR,'L')
            handles.RectL = rectangle('Position',handles.tempRect1);
            set(handles.RectL,'EdgeColor','r','LineWidth',3);
        elseif strcmp(handles.LorR,'R')
            handles.RectR = rectangle('Position',handles.tempRect1);
            set(handles.RectR,'EdgeColor','r','LineWidth',3);
        end
    end
    guidata(hObject,handles);
end



% --- Executes on button press in Reset.
function Reset_Callback(hObject, eventdata, handles)
% hObject    handle to Reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
handles.checkbox1.Value=0;
cla;
imagesc(handles.mI);
colormap gray;
hold on;
plot(handles.x,handles.y,'k.-');

handles.RectL = [];
handles.RectR = [];
handles.Line = [];
guidata(hObject,handles);



function AreaVal_Callback(hObject, eventdata, handles)
% hObject    handle to AreaVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AreaVal as text
%        str2double(get(hObject,'String')) returns contents of AreaVal as a double


% --- Executes during object creation, after setting all properties.
function AreaVal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AreaVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in EdgeBox.
function EdgeBox_Callback(hObject, eventdata, handles)
% hObject    handle to EdgeBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of EdgeBox
handles.P1 = [];
handles.P2 = [];
handles.P3 = [];
handles.P4 = [];
handles = guidata(hObject);

if get(hObject,'Value')==1
    
    handles.EdgeThrVal = str2double(get(handles.EdgeThr,'String'));
    
    pL = handles.RectL.Position;
    if ~isempty(handles.RectR)
        pR = handles.RectR.Position;
    else
        pR = [];
    end
    pT = handles.EdgeThrVal*10*handles.pixpermm;
    x = handles.x;
    y = handles.y;
    
    if ~isempty(handles.Line)
        indxL = x < mean(handles.Line.XData);
        indxR = x > mean(handles.Line.XData);
    else
        indxL = ones(size(x));
        indxR = zeros(size(x));
        warning('assuming all position data is in LEFT box');
    end
    
    edgeL = (x<pL(1)+pT | y<pL(2)+pT | x>pL(1)+pL(3)-pT | y>pL(2)+pL(4)-pT) & indxL;
    if ~isempty(pR)
        edgeR = (x<pR(1)+pT | y<pR(2)+pT | x>pR(1)+pL(3)-pT | y>pR(2)+pL(4)-pT) & indxR;
    else
        edgeR = dec2bin(zeros(size(x)));
    end
    centL = ~edgeL & ~indxR;
    centR = ~edgeR & ~indxL;
    
    handles.P1 = plot(handles.x(edgeL),handles.y(edgeL),'b.');
    handles.P2 = plot(handles.x(centL),handles.y(centL),'c.');
    handles.P3 = plot(handles.x(edgeR),handles.y(edgeR),'b.');
    handles.P4 = plot(handles.x(centR),handles.y(centR),'c.');
    
    handles.edgeL = edgeL;
    handles.edgeR = edgeR;
    handles.centL = centL;
    handles.centR = centR;
    
    guidata(hObject,handles);
else
    delete(handles.P1);
    delete(handles.P2);
    delete(handles.P3);
    delete(handles.P4);
end



% --- Executes on button press in Exclusion_Zone.
function Exclusion_Zone_Callback(hObject, eventdata, handles)
% hObject    handle to Exclusion_Zone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.ExclusionZone = getrect;
handles.ExclusionZone = rectangle('Position',handles.ExclusionZone);
set(handles.ExclusionZone,'EdgeColor','g','LineWidth',2);

guidata(hObject,handles);



function EdgeThr_Callback(hObject, eventdata, handles)
% hObject    handle to EdgeThr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EdgeThr as text
%        str2double(get(hObject,'String')) returns contents of EdgeThr as a double

handles.EdgeThrVal = str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function EdgeThr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EdgeThr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in crossings.
function crossings_Callback(hObject, eventdata, handles)
% hObject    handle to crossings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of crossings

if get(hObject,'Value')==1
    
    delete(handles.P1);
    delete(handles.P2);
    delete(handles.P3);
    delete(handles.P4);
    
    pL = handles.RectL.Position;
    
    if ~isempty(handles.RectR)
        
        delete(handles.P1);
        delete(handles.P2);
        delete(handles.P3);
        delete(handles.P4);
        
        pR = handles.RectR.Position;
        pT = handles.EdgeThrVal*10*handles.pixpermm;
        x = handles.x;
        y = handles.y;
        
        % compute crossingls
        inL = x < pL(1)+pL(3);
        inR = x > pR(1);
        
        %keyboard;
        
        iR = find(cat(1,0,abs(diff(inR))>0));
        iL = find(cat(1,0,abs(diff(inL))>0));
        ixR = ones(size(iR));
        ixL = ones(size(iL))+2;
        X = cat(1,ixR,ixL);
        [S,idx] = sort(cat(1,iR,iL));
        iX = X(idx);
        
        nX = length(idx);
        newix = iX;
        newS = S;
        cix = [];
        nix = [];
        % get first ix
        firstix = newix(1);
        switch firstix
            case 1 % R
                currix = 1;
                nextix = 3;
            case 3 % L
                currix = 3;
                nextix = 1;
        end
        lastchange = find(diff(newix)>0,1,'last');
        endS = newS(lastchange+1);
        switch newix(lastchange+1)
            case 1
                lastix = 1;
            case 3
                lastix = 3;
        end
        cix = cat(1,cix,currix);
        nix = cat(1,nix,nextix);
        for ii=2:nX-1
            temp = newix(ii:end);
            if currix==nextix
                switch currix
                    case 1
                        nextix = 3;
                    case 3
                        nextix = 1;
                end
            elseif currix~=nextix
                currix = nextix;
            end
            cix = cat(1,cix,currix);
            nix = cat(1,nix,nextix);
            
            
            if currix==1 && nextix==1 || currix==3 && nextix==1
                temp = find(temp==1,1,'first');
            elseif currix==1 && nextix==3 || currix==3 && nextix==3
                temp = find(temp==3,1,'first');
            end
            if ~isempty(temp)
                newix(ii:temp)=[];
                newS(ii:temp)=[];
            else
                newix(ii:end)=[];
                newS(ii:end)=[];
            end
            if ii+1>length(newix)
                break;
            end
        end
        
        % find 'repeats' of length > 3 and get 'edges'
        
        nx = length(newix);
        sumR = 1;
        while sumR>0
            sumR = 0;
            for ix=2:nx-1
                if newix(ix)==newix(ix-1) && newix(ix)==newix(ix+1)
                    newS(ix)=[];
                    newix(ix)=[];
                    sumR = sumR+1;
                end
                if ix+1 > length(newix)-1
                    break;
                end
            end
        end
        
        if newix(end)==lastix
        else
            newix(end+1) = lastix;
            newS(end+1) = endS;
        end
        
        handles.nRL = sum(diff(newix)==2);
        handles.nLR = sum(diff(newix)==-2);
        handles.LtoR = zeros(size(inR));
        handles.RtoL = zeros(size(inR));
        handles.LtoL = zeros(size(inR));
        handles.RtoR = zeros(size(inR));
        handles.tRL = nan(handles.nRL,1);
        handles.tLR = nan(handles.nLR,1);
        nx = length(newix);
        
        switch firstix
            case 1
                handles.RtoR(1:newS(1))=1;
            case 3
                handles.LtoL(1:newS(1))=1;
        end
        iC = 1;
        iD = 1;
        for ix=2:nx
            switch newix(ix-1)
                case 1
                    switch newix(ix)
                        case 1 % R->R
                            handles.RtoR(newS(ix-1):newS(ix))=1;
                        case 3 % R->L
                            handles.RtoL(newS(ix-1):newS(ix))=1;
                            handles.tRL(iC,1)=newS(ix-1)./30;
                            iC = iC+1;
                    end
                case 3
                    switch newix(ix)
                        case 1 % L->R
                            handles.LtoR(newS(ix-1):newS(ix))=1;
                            handles.tLR(iD,1)=newS(ix-1)./30;
                            iD = iD+1;
                        case 3 % L->L
                            handles.LtoL(newS(ix-1):newS(ix))=1;
                    end
            end
        end
        switch lastix
            case 1
                handles.RtoR(newS(ix):end)=1;
            case 3
                handles.LtoL(newS(ix):end)=1;
        end
        delete(handles.P1);
        delete(handles.P2);
        delete(handles.P3);
        delete(handles.P4);
        
        handles.inL = inL;
        handles.inR = inR;
        
        k = ~handles.LtoR & ~handles.RtoL;
        handles.P1 = plot(x(handles.edgeL & k),y(handles.edgeL & k),'c.');
        handles.P2 = plot(x(handles.edgeR & k),y(handles.edgeR & k),'c.');
        handles.P3 = plot(x(handles.centL & k),y(handles.centL & k),'b.');
        handles.P4 = plot(x(handles.centR & k),y(handles.centR & k),'b.');
        
        handles.binSize = str2double(get(handles.binVal,'String'));
        
        nP = length(inL);
        nB = linspace(1,nP,ceil(nP./(30*handles.binSize)));
        tL_binned = histcounts(find(handles.inL),nB);
        tR_binned = histcounts(find(handles.inR),nB);
        edgeL_binned = histcounts(find(handles.edgeL),nB);
        edgeR_binned = histcounts(find(handles.edgeR),nB);
        centL_binned = histcounts(find(handles.centL),nB);
        centR_binned = histcounts(find(handles.centR),nB);
        handles.tL_binned = tL_binned/30;
        handles.tR_binned = tR_binned/30;
        handles.edgeL_binned = edgeL_binned/30;
        handles.edgeR_binned = edgeR_binned/30;
        handles.centL_binned = centL_binned/30;
        handles.centR_binned = centR_binned/30;
        
        guidata(hObject,handles);
    else
        
        x = handles.x;
        y = handles.y;
        
        % compute crossingls
        inL = x < pL(1)+pL(3);
        inR = zeros(size(x));
        handles.inL = inL;
        handles.inR = inR;
        
        %k = ~handles.LtoR & ~handles.RtoL;
        handles.P1 = plot(x(handles.edgeL),y(handles.edgeL),'c.');
        handles.P2 = plot(x(handles.edgeR),y(handles.edgeR),'c.');
        handles.P3 = plot(x(handles.centL),y(handles.centL),'b.');
        handles.P4 = plot(x(handles.centR),y(handles.centR),'b.');
        
        handles.binSize = str2double(get(handles.binVal,'String'));
        
        nP = length(inL);
        nB = linspace(1,nP,ceil(nP./(30*handles.binSize)));
        tL_binned = histcounts(find(handles.inL),nB);
        tR_binned = histcounts(find(handles.inR),nB);
        edgeL_binned = histcounts(find(handles.edgeL),nB);
        edgeR_binned = histcounts(find(handles.edgeR),nB);
        centL_binned = histcounts(find(handles.centL),nB);
        centR_binned = histcounts(find(handles.centR),nB);
        handles.tL_binned = tL_binned/30;
        handles.tR_binned = tR_binned/30;
        handles.edgeL_binned = edgeL_binned/30;
        handles.edgeR_binned = edgeR_binned/30;
        handles.centL_binned = centL_binned/30;
        handles.centR_binned = centR_binned/30;
        
        guidata(hObject,handles);
    end
else
    delete(handles.P1);
    delete(handles.P2);
    delete(handles.P3);
    delete(handles.P4);
end




function binVal_Callback(hObject, eventdata, handles)
% hObject    handle to binVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of binVal as text
%        str2double(get(hObject,'String')) returns contents of binVal as a double


% --- Executes during object creation, after setting all properties.
function binVal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to binVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
