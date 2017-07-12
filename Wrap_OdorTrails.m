function [x,y,V,dT,Xp,Yp,xT1,yT1,idphi,nidphi,mouse,trial,sess,conc,frame,zidphi,znidphi,nearX,nearY,nx,ny,nV,mouseT,sessT,mouse1,sess1,mouseName1,mouseName2,dnT,Xnp,Ynp,C,nC,znC,zC] = Wrap_OdorTrails
%function [x,y,V,dT,Xp,Yp,xT1,yT1,idphi,nidphi,mouse,trial,sess,conc,frame,zidphi,znidphi,nearX,nearY,nx,ny,nV,mouseT,sessT,mouse1,sess1,mouseName2,C,nC,znC,zC] = Wrap_OdorTrails
% 2017-07-03 AndyP
% Wrapper function to process and concatenate position data from optimouse
% and extracted trails from getTrail_GUI3.  Position data and derivatives are in the format (nframes x
% 1) where nframes is the total number of frames to be analyzed.  Trail
% data is in the format (npoints x 1) where n points is the total number of
% points for all trails.
% CONVERSIONS
% pixels to cm = 11.2 pixels / cm
% INPUTS
% readXLS 1x1 logical.
% If true, mouse, sess, conc, trial, mouseT, and sessT are taken from an xls
% worksheet with notes in columns 3-6 about the type of trail.  The mouse
% names are in column 2 and the dates are in column 1.  Assumes dates are
% all in a standardized format with dates YYDDMM. (Y=year, D=day, M=month).
% If false, conc and trial are empty and mouse, sess, mouseT, and sessT are taken from parsing the filenames in the for loop.
% mouse1 and sess1 are checks on mouse and sess. These assertions should be
% true: all(mouse==mouse1); all(sess==sess1);
% pathType  string
% 'Y' 'curvy' 'zig zag' or 'spot'
% OUTPUTS
% x  - (nframes x 1 double) <x> positions from optimouse processed with
% Process_VT
% y  - (nframes x 1 double) <y> positions from optimouse processed with
% Process_VT
% V - (nframes x 1 double) [cm/s] velocity sqrt(dy.^2+dx.^2) using dxdt function (Janabi-Sharifi
% algorithm).
% dT - (nframes x 1 double) [cm] distance from current <x,y> position to closest
% trail coordilate
% xP - (nframes x 1 double) x position trail coordinate closest to <x,y>
% yP - (nframes x 1 double) y position trail coordinate closest to <x,y>
% xT1 - (npoints x 1 double) x position trail coordinates
% yT1 - (npoints x 1 double) y position trail coordinates
% idphi - derivative of tangent angle between <dx,dy> computed using dxdt
% function (Janabi-Sharifi algorithm).
% nidphi - derivative of tangent angle of nose, computed usuing dxdt
% function (Janabi-Sharifi algorithm).
% mouse - (nframes x 1 double) mouse counter (starts at 1, +1 for each new mouse)
% trial - (nframes x 1 double) trial counter (starts at 1, +1 for each
% new trial, typically 4-5 trials/session)
% sess - (nframes x 1 double) session counter (starts at 1, +1 for each new
% session, typically 1 session/day with 4-5 trials/session).
% conc - (nframes x 1 double) odor concentration for the given session,
% derived from the xls worksheet notes.  [] if readXLS is false.
% frame - (nframes x 1 double) the current frame of the position sample
% (typically 9000 frames / trial at 50 frames / s).
% zidphi - (nframes x 1 double) z-scored absolute value of idphi
% znidphi - (nframes x 1 double) z-scored absolute value of nidphi
% nearX - (nframes x 1 double) interpolated <x> positions from optimouse, using nearest non-nan
% x-position.  sum(isnan(nearX)) should be true.
% nearY - (nframes x 1 double) interpolated <y> positions from optimouse, using nearest non-nan
% x-position.  sum(isnan(nearX)) should be true.
% nx - (nframes x 1 double) nose <x> position from optimouse
% ny - (nframes x 1 double) nose <y> position from optimouse
% nV - (nframes x 1 double) [cm/s] nose velocity (sqrt(dny.^2+dnx.^2))
% using dxdt function (Janabi-SHarifi algorithm).
% mouseT - (npoints x 1 double) mouse counter for the trail points <xT1,yT1> (starts
% at 1, +1 for each new mouse).
% sessT - (npoints x 1 double) session counter for the trail points
% <xT1,yT1> (starts at 1, +1 for each new session).
% mouse1 - should be identical to mouse
% sess1 - should be identical to sess
% mouseName1 - (nS x 1 double) mouse name counter
% mouseName2 - should be identical to mouseName1
% dnT - (nframes x 1 double) [cm] nose distance from current <nx,ny>
% position to closest trail coordinate
% Xnp - (nframes x 1 double) closest trail <x> coordinate to nose
% Ynp (nframes x 1 double) closest trail <y> coordinate to nose
% EXAMPLE
%After running the wrapper, to extract position samples from a single sesion,
% use the mouse and sess outputs.
% nM = length(unique(mouse)); nS=max(sess);
% for iM=1:nM
%     for iS=1:nS
%         currSess = mouse==iM & sess==iS;
%         currTrail = mouseT==iM & sessT==iS;
%         plot(xT(currTrail),yT(currTrail),'r.'); % plot current trail
%         hold on;
%         plot(x(currSess),y(currSess),'.','markersize',1,'color',[0.5 0.5 0.5]); % plot position samples for current session
%         pause;
%     end
% end


homedir = cd;

readXLS = true;
pathType = 'Y';
postSmoothing = 0.1; % s
window = 1; % s
dtime = 1/50; % Hz

trial0 = [];
sess0 = [];
conc0 = [];
mouse0 = [];
mouseName1 = [];
iS0 = 1;
iM = 0;
lastMouse = [];
dateStr = [];

if readXLS
    
    % get trial numbers from xls sheet...
    warning('need to change directory and file name to match computer directory location for xls notesheet');
    cd('C:\Users\papalea\Documents\Data');
    [num,text] = xlsread('Log_170427.xlsx',1,'A1:H176','basic'); % read xls sheet into matlab
    
    for iS=2:size(text,1);
        mousetemp = num(iS-1,2); % mice are always in second column, starting from row 1
        
        for iT=4:8
            tempStr = text{iS,iT};
            switch pathType
                case 'Y'
                    k = strfind(tempStr,'y');
                    k1 = strfind(tempStr,'Y');
                    k4 = strfind(tempStr,'lowercase y');
                    k2 = strfind(tempStr,'curvy');
                    k3 = strfind(tempStr,'Curvy');
                case 'curvy'
                    k4 = strfind(tempStr,'curvy'); % k4 anywhere
                    k1 = strfind(tempStr,'Curvy'); % k1<5
                    k = strfind(tempStr,'curved'); % k<5
                    k2 = []; % if for some reason there is curvy,Curvy,curved in non-curvy trials
                    k3 = []; % ditto to k2
                case 'zig zag'
                    k4 = strfind(tempStr,'zig zag'); % k4 anywhere
                    k1 = []; % can also be true, k1 < 5
                    k = [];  % can also be true, k < 5
                    k2 = []; % false positive, if for some reason there is zig zag in non zig zag trials
                    k3 = []; % false positive, see k2
                case 'spot'
                    k4 = strfind(tempStr,'spot'); % see above comments
                    k1 = [];
                    k = [];
                    k2 = [];
                    k3 = [];
                otherwise
                    error('unknown pathType');
            end
            if (~isempty(k) & k<5 | ~isempty(k1) & k1<5 | ~isempty(k4)) & isempty(k2) & isempty(k3) %#ok<OR2,AND2>
                mouseName1 = cat(1,mouseName1,mousetemp);
                dateStr = cat(1,dateStr,num(iS-1,1));
                %disp(mouse0);
                disp(tempStr);
                trial0 = cat(1,trial0,iT-3); % trial 1 is in column 4, trial 2 in column 5, etc.
                % code to update mouse0 / sess0
                if lastMouse==mousetemp
                    iS0=iS0+1;
                else
                    iS0=1;
                    iM=iM+1;
                end
                mouse0 = cat(1,mouse0,iM);
                sess0 = cat(1,sess0,iS0);
                lastMouse = mousetemp;
                % code to get concentration from xls text
                kone = strfind(tempStr,'1%');
                ktwo = strfind(tempStr,'2%');
                if ~isempty(kone) && isempty(ktwo)
                    conc0 = cat(1,conc0,1);
                elseif isempty(kone) && ~isempty(ktwo)
                    conc0 = cat(1,conc0,2);
                else
                    warning('unknown concentration');
                    conc0 = cat(1,conc0,nan);
                end
            end
        end
    end
    
    cd(homedir);
    
end

% get position and trail files
positFiles = dir('*_positions.mat');
trailFiles = dir('*-Odor-Trail.mat');
startFiles = dir('*-StartFrame.mat');

assert(length(positFiles)==length(trailFiles),'error needs to be an odor trail file for each position file');
assert(length(startFiles)==length(trailFiles),'error needs to be a start time for each position file');

if readXLS
    
    % find missing positFiles
    for iT=1:length(mouseName1);
        currName = mouseName1(iT);
        currDate = dateStr(iT);
        flag = 0;
        for iF=1:length(positFiles)
            tempStr = strsplit(positFiles(iF).name,'_');
            mouseStr = tempStr{1};
            dateStr1 = tempStr{2}(1:10);
            yr = dateStr1(3:4);
            mo = dateStr1(6:7);
            da = dateStr1(9:10);
            dateStr1 = strcat(yr,mo,da);
            if isequal(currDate,str2double(dateStr1)) && isequal(currName,str2double(mouseStr))
                flag = 1;
            else
            end
        end
        if flag
        else
            disp(currName); disp(currDate);
        end
    end
    
    % find missing positFiles
    for iT=1:length(positFiles);
        tempStr = strsplit(positFiles(iT).name,'_');
        mouseStr = tempStr{1};
        dateStr1 = tempStr{2}(1:10);
        yr = dateStr1(3:4);
        mo = dateStr1(6:7);
        da = dateStr1(9:10);
        dateStr1 = strcat(yr,mo,da);
        
        flag = 0;
        for iF=1:length(mouseName1)
            currName = mouseName1(iF);
            currDate = dateStr(iF);
            if isequal(currDate,str2double(dateStr1)) && isequal(currName,str2double(mouseStr))
                flag = 1;
            else
            end
        end
        if flag
        else
            disp(tempStr);
        end
    end
    
    
    assert(length(positFiles)==length(trial0),'error notes in xls sheet missing entry or missing position file');
    
end

nD = length(trailFiles);

x = [];
y = [];
V = [];
xT1 = [];
yT1 = [];
dT = [];
Xp = [];
Yp = [];
C = [];
nC = [];
zC = [];
znC = [];
idphi = [];
nidphi = [];
mouse = [];
trial = [];
sess = [];
conc = [];
frame = [];
nearX = [];
nearY = [];
nx = [];
ny = [];
nV = [];
mouseT = [];
sessT = [];
mouse1 = [];
sess1 = [];
lastMouse = '';
iM = 1;
iS = 0;
mouseName2 = [];
zidphi = [];
znidphi = [];
dnT = [];
Xnp = [];
Ynp = [];
for iD=1:nD
    cd(homedir);
    fprintf('%s %d/%d \n',trailFiles(iD).name,iD,nD);
    load(positFiles(iD).name,'position_results');
    load(trailFiles(iD).name,'xT','yT','data');
    load(startFiles(iD).name,'startFrame');
    
    % code to process the position output from optimouse
    [x0,y0,nx0,ny0] = Process_VT(position_results,startFrame);
    
    % code to get mouse1 / sess1, should be identical to mouse/sess from
    % xls file
    tempStr = strsplit(positFiles(iD).name,'_');
    mouseStr = tempStr{1};
    mouseName2 = cat(1,mouseName2,str2double(mouseStr));
    if strcmp(mouseStr,lastMouse) | iD==1 %#ok<OR2>
        iS=iS+1;
    else
        iM=iM+1;
        iS=1;
    end
    lastMouse = mouseStr;
    mouse1 = cat(1,mouse1,repmat(iM,[length(x0),1]));
    sess1 = cat(1,sess1,repmat(iS,[length(x0),1]));
    
    
    % Get trail coordinates.  Note y,x correctly switched from trailFile data.
    xT0 = yT;
    yT0 = xT;
    
    xT1 = cat(1,xT1,xT0);
    yT1 = cat(1,yT1,yT0);
    
    
    % for each nan, find nearest position to each nan...
    
    %// Index array for factor
    x1 = 1:numel(x0);
    %// Indices of NaNs
    t2 = find(~isnan(x0));
    
    %// Replace NaNs with the closest non-NaNs
    nearX0 = interp1(x1(t2),x0(t2),x1,'nearest');
    nearY0 = interp1(x1(t2),y0(t2),x1,'nearest');
    
    nearX0(1)=nearX0(2);
    nearY0(1)=nearY0(2);
    
    nearX = cat(1,nearX,nearX0');
    nearY = cat(1,nearY,nearY0');
    
    % get position data
    x = cat(1,x,x0);
    y = cat(1,y,y0);
    nx = cat(1,nx,nx0);
    ny = cat(1,ny,ny0);
    
    % get frame data
    frame = cat(1,frame,(1:length(x0))');
    
    % compute velocity
    dx = dxdt(x0,dtime,window,postSmoothing);
    dy = dxdt(y0,dtime,window,postSmoothing);
    V0 = sqrt(dx.^2+dy.^2)./11.2;
    V = cat(1,V,V0); % cm/s
    
    ndx = dxdt(nx0,dtime,window,postSmoothing);
    ndy = dxdt(ny0,dtime,window,postSmoothing);
    nV0 = sqrt(ndx.^2+ndy.^2)./11.2;
    nV = cat(1,nV,nV0); % cm/s
    
    % calculate distance from trail
    nP = length(x0);
    dT0 = nan(nP,1);
    I = nan(nP,1);
    dnT0 = nan(nP,1);
    In = nan(nP,1);
    for iP=1:nP
        [dT0(iP),I(iP)] = nanmin(sqrt((x0(iP)-xT0).^2+(y0(iP)-yT0).^2));
        [dnT0(iP),In(iP)] = nanmin(sqrt((nx0(iP)-xT0).^2+(ny0(iP)-yT0).^2));
    end
    dT = cat(1,dT,dT0./11.2);
    Xp = cat(1,Xp,xT0(I));
    Yp = cat(1,Yp,yT0(I));
    dnT = cat(1,dnT,dnT0./11.2);
    Xnp = cat(1,Xnp,xT0(In));
    Ynp = cat(1,Ynp,yT0(In));
    
    
    C0 = Tortuosity(x0,y0,dtime,window,postSmoothing);
    C = cat(1,C,C0);
    %
    nC0 = Tortuosity(nx0,ny0,dtime,window,postSmoothing);
    nC = cat(1,nC,nC0);
    
    % compute idphi
    idphi0 = zIdPhi(dx,dy,dtime,window,postSmoothing);
    nidphi0 = zIdPhi(ndx,ndy,dtime,window,postSmoothing);
    
    idphi = cat(1,idphi,idphi0);
    nidphi = cat(1,nidphi,nidphi0);
    
    zidphi= cat(1,zidphi,nanzscore(abs(idphi0)));
    znidphi= cat(1,znidphi,nanzscore(abs(nidphi0)));
    
    zC = cat(1,zC,nanzscore(abs(C0)));
    znC = cat(1,znC,nanzscore(abs(nC0)));
    
    if readXLS
        mouse = cat(1,mouse,repmat(mouse0(iD),[length(x0),1])); %#ok<UNRCH>
        trial = cat(1,trial,repmat(trial0(iD),[length(x0),1]));
        conc = cat(1,conc,repmat(conc0(iD),[length(x0),1]));
        sess = cat(1,sess,repmat(sess0(iD),[length(x0),1]));
        mouseT = cat(1,mouseT,repmat(mouse0(iD),[length(xT0),1]));
        sessT = cat(1,sessT,repmat(sess0(iD),[length(xT0),1]));
    else
        mouse = cat(1,mouse,repmat(iM,[length(x0),1]));
        sess = cat(1,sess,repmat(iS,[length(x0),1]));
        mouseT = cat(1,mouseT,repmat(iM,[length(xT0),1]));
        sessT = cat(1,sessT,repmat(iS,[length(xT0),1]));
    end
    
end

end

