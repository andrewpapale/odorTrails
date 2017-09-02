function [x,y,V,dT,Xp,Yp,idphi,nidphi,mouse,trial,sess,conc,frame,zidphi,znidphi,nearX,nearY,nx,ny,nV] = Wrap_Generic
%[x,y,V,dT,Xp,Yp,idphi,nidphi,mouse,trial,sess,conc,frame,zidphi,znidphi,nearX,nearY,nx,ny,nV] = Wrap_Generic
% 2017-07-03 AndyP
% 2017-07-26 AndyP updated with foaw_diff, Tortuosity1 and zIdPhi1
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

postSmoothing = 0.25; % s
%window = 0.5; % s
m = 50;
d = 1;
dtime = 1/30; % Hz


% get position and trail files
positFiles = dir('*_positions.mat');
startFiles = dir('*-StartFrame.mat');

nD = length(positFiles);

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
bait = [];
for iD=1:nD
    cd(homedir);
    fprintf('%s %d/%d \n',positFiles(iD).name,iD,nD);
    load(positFiles(iD).name,'position_results');
    load(startFiles(iD).name,'startFrame');
    
    %    code to process the position output from optimouse
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
    dx = foaw_diff(x0, dtime, m, d, postSmoothing);
    dy = foaw_diff(y0, dtime, m, d, postSmoothing);
    V0 = sqrt(dx.^2+dy.^2)./11.2;
    V = cat(1,V,V0); % cm/s
    
    ndx = foaw_diff(nx0, dtime, m, d, postSmoothing);
    ndy = foaw_diff(ny0, dtime, m, d, postSmoothing);
    nV0 = sqrt(ndx.^2+ndy.^2)./11.2;
    nV = cat(1,nV,nV0); % cm/s
    
    C0 = Tortuosity1(dx,dy,dtime,m,d,postSmoothing);
    C = cat(1,C,C0);
    %
    nC0 = Tortuosity1(ndx,ndy,dtime,m,d,postSmoothing);
    nC = cat(1,nC,nC0);
    
    % compute idphi
    idphi0 = zIdPhi1(dx,dy,dtime,m,d,postSmoothing);
    nidphi0 = zIdPhi1(ndx,ndy,dtime,m,d,postSmoothing);
    
    idphi = cat(1,idphi,idphi0);
    nidphi = cat(1,nidphi,nidphi0);
    
    zidphi= cat(1,zidphi,nanzscore(abs(idphi0)));
    znidphi= cat(1,znidphi,nanzscore(abs(nidphi0)));
    
    zC = cat(1,zC,nanzscore(abs(C0)));
    znC = cat(1,znC,nanzscore(abs(nC0)));
    
    
    
    mouse = cat(1,mouse,repmat(iM,[length(x0),1]));
    sess = cat(1,sess,repmat(iS,[length(x0),1]));
end
end

