function [usedLEAP,x,y,nx,ny,tx,ty,tbx,tby,confb,confn,conftb,conftt,cx1,cy1,t,mL,sess,mouse,dT,dnT,nTr,t2f,t2s,f2s,f2f,sessT,mouseT,sT,t2s1,mTlen,mTlen2,C,sP,Tr,sessS,mouseS,initD,orient,initO,bV,nV,ndphi,bdphi,slope,Dline,theta,edge,cx0,cy0] = Wrap_autoSpot()


homedir = cd;
readXLS = true;

% 1 x (nP*nT*nS)
usedLEAP = [];
x = [];
y = [];
nx1 = [];
ny1 = [];
tx = [];
ty = [];
tbx1 = [];
tby1 = [];
cx1 = [];
cy1 = [];
confb1 = [];
confn1 = [];
conftb1 = [];
conftt1 = [];
t = [];
mL = [];
orient = [];
bV = [];
nV = [];
ndphi = [];
bdphi = [];
slope = [];
Dline = [];
theta = [];

sess = [];
mouse = [];
trial = [];
dT = [];
dnT = [];

% 1x (nT*nS)
nTr = [];
t2f = [];
t2s = [];
f2s = [];
f2f = [];
mouseT = [];
sessT = [];

% get mice in directory
D = dir;
nD = length(D);
nM = 0;
iC = 1;
for iD=1:nD
    if D(iD).isdir && ~strcmp(D(iD).name,'.') && ~strcmp(D(iD).name,'..') && ~strcmp(D(iD).name,'models') && ~strcmp(D(iD).name,'old data')
        nM = nM+1;
        mDir{iC} = D(iD).name; %#ok<AGROW>
        iC=iC+1;
    end
end

for iM=1:nM
    cd(mDir{iM});
    % get sessions
    D = dir('*.mat');
    nD = length(D); % includes LEAP and matlab-generated mat files
    ismat = [];
    for iD=1:nD
        matexp = '\d\d\d\d-\d\d-\d\d_\d\d_\d\d.mat';
        leapexp = '\d\d\d\d\d\d-\d\d\d\d-\d\d-\d\d-\d.mat';
        ismat0 = regexp(D(iD).name,matexp,'once');
        if ~isempty(ismat0) % leap file
            ismat = cat(1,ismat,ismat0);
        else
            ismat = cat(1,ismat,0);
        end
    end
    
    D1 = D;
    L1 = D;
    D1(~ismat)=[];
    L1(ismat==1)=[];
    
    nD = length(D1);
    for iD=1:nD
        sList{1,iM}{iD} = D1(iD).name;
        sList{2,iM}{iD} = D1(iD).name(1:10);
    end
    
    nL = length(L1);
    for iL=1:nL
        sListL{1,iM}{iL} = L1(iL).name;
        sListL{2,iM}{iL} = L1(iL).name(8:17);
    end
    
    cd(homedir);
end

nS = 105;
% (nM x nS)
mouseS = nan(nM,nS);
sessS = nan(nM,nS);
sT = nan(nM,nS);  % total number of trials
t2s1 = nan(nM,nS); % first time to spot (s)
mTlen = nan(nM,nS); % mean trial length
mTlen2= nan(nM,nS); % median trial length
C = nan(nM,nS); % concentration
sP = nan(nM,nS); % spot position (A-N -> 1-14)
Tr = nan(nM,nS); % trial (1-5)
initD = nan(nM,nS);
initO = nan(nM,nS);


if readXLS
    [num,text] = xlsread('ConcentrationExperiment_2018-10-22.xlsx','A1:E89');
    
    % for each day iD, get mouse num(iD,1) and date num(iD,2)
    % get number of sessions per day length(strsplit(text{iD,2},','))
    % get number of .mat files for that mouse and day
    % make sure number of .mat files = number of sessions/day, else error
    % read text data spot ID, % ID, experimenter ID
    
    nD = size(num,1);
    for iD=1:nD
        miD = mat2str(num(iD,1));
        % switch mouse
        if iD==1
            lastmiD = 'dum';
        else
            lastmiD = mat2str(num(iD-1,1));
        end
        if ~strcmp(miD,lastmiD)
            iC = 1;
        else
        end
        
        
        diD = mat2str(num(iD,2));
        piD = strsplit(text{iD,3},',');
        siD = strsplit(text{iD,2},',');
        
        nansiD = zeros(size(siD));
        for iCell=1:length(siD)
            nansiD(iCell)=strcmp(siD{iCell},'nan');
        end
        
        assert(length(siD)==length(piD),'xls error: mismatch in Spot ID and % ID numbers');
        tempA = cellfun(@strcmp,mDir,repmat({miD},[1,nM]));
        assert(sum(tempA)==1,'error unknown mouse, likely xls typo');
        
        ystr = diD(1:2);
        mstr = diD(3:4);
        dstr = diD(5:6);
        dateStr = strcat('20',ystr,'-',mstr,'-',dstr);
        
        mI = find(tempA==1); % mouse Index
        nS0 = length(sList{1,mI});
        tempB = cellfun(@strcmp,sList{2,mI},repmat({dateStr},[1,nS0]));
        assert(sum(tempB)==length(siD)-sum(nansiD),'xls/mat file mismatch for %d, %s',iD,dateStr);
        
        % check that for each D1 there is an L1, otherwise fill with a
        % blank
        
        nSL = length(sListL{1,mI});
        tempL = cellfun(@strcmp,sListL{2,mI},repmat({dateStr},[1,nSL]));
        
        % for debugging
        %         disp(sList{1,mI}(tempB==1)');
        %         disp(sListL{1,mI}(tempL==1)');
        
        sListL0 = [];
        if sum(tempL)~=sum(tempB)
            % if vid is missing, determine which video is missing
            % if extra vid, determine which vid is extra
            L = sum(tempL);
            S = sum(tempB);
            
            Lsess = find(tempL==1);
            
            for iL=1:length(Lsess)
                sListL0{iL,1} = sListL{1,mI}{Lsess(iL)};
                sListL0{iL,2} = Lsess(iL);
                sListL0{iL,3} = sListL{1,mI}{Lsess(iL)}(19);
            end
            
            sessNum = cellfun(@str2double,sListL0(:,3))';
            
            
            if L > S % extra video, must be at least one removed session
                
            elseif S > L % vid missing, must add at least 1 space to sListL
                if length(sessNum)==4 && sessNum(1)==1 && sessNum(2)==2 && sessNum(3)==3 && sessNum(4)==4 && ~strcmp(dateStr,'2018-08-30')
                    sListL{1,mI} = [sListL{1,mI}(1:Lsess(4)),nan,sListL{1,mI}(Lsess(4)+1:end)];
                    sListL{2,mI} = [sListL{2,mI}(1:Lsess(4)),nan,sListL{2,mI}(Lsess(4)+1:end)];
                elseif length(sessNum)==4 && sessNum(1)==1 && sessNum(2)==2 && sessNum(3)==3 && sessNum(4)==5
                    sListL{1,mI} = [sListL{1,mI}(1:Lsess(3)),nan,sListL{1,mI}(Lsess(4):end)];
                    sListL{2,mI} = [sListL{2,mI}(1:Lsess(3)),nan,sListL{2,mI}(Lsess(4):end)];
                elseif length(sessNum)==4 && sessNum(1)~=1
                    sListL{1,mI} = [sListL{1,mI}(1:Lsess(1)-1),nan,sListL{1,mI}(Lsess(1):end)];
                    sListL{2,mI} = [sListL{2,mI}(1:Lsess(1)-1),nan,sListL{2,mI}(Lsess(1):end)];
                elseif length(sessNum)==3 && sessNum(1)==1 && sessNum(2)==2 && sessNum(3)==3
                    sListL{1,mI} = [sListL{1,mI}(1:Lsess(3)),nan,nan,sListL{1,mI}(Lsess(3)+1:end)];
                    sListL{2,mI} = [sListL{2,mI}(1:Lsess(3)),nan,nan,sListL{2,mI}(Lsess(3)+1:end)];
                elseif length(sessNum)==3 && sessNum(1)~=1 && sessNum(2)==4 && strcmp(dateStr,'2018-07-24')
                    sListL{1,mI} = [nan,sListL{1,mI}(Lsess(1):end)];
                    sListL{2,mI} = [nan,sListL{2,mI}(Lsess(1):end)];
                elseif length(sessNum)==4 && sessNum(1)==1 && sessNum(2)==3
                    sListL{1,mI} = [sListL{1,mI}(1:Lsess(1)),nan,sListL{1,mI}(Lsess(2):end)];
                    sListL{2,mI} = [sListL{2,mI}(1:Lsess(1)),nan,sListL{2,mI}(Lsess(2):end)];
                elseif length(sessNum)==3 && sessNum(1)==2 && sessNum(2)==4 && strcmp(dateStr,'2018-07-10')
                    sListL{1,mI} = [nan,sListL{1,mI}(Lsess(1)),nan,sListL{1,mI}(Lsess(2)),sListL{1,mI}(Lsess(3):end)];
                    sListL{2,mI} = [nan,sListL{2,mI}(Lsess(1)),nan,sListL{2,mI}(Lsess(2)),sListL{2,mI}(Lsess(3):end)];
                elseif length(sessNum)==4 && sessNum(1)==1 && sessNum(2)==2 && sessNum(3)==3 && sessNum(4)==4 && strcmp(dateStr,'2018-08-30')
                    sListL{1,mI} = [sListL{1,mI}(1:Lsess(4)),nan,sListL{1,mI}(Lsess(4)+1:end)];
                    sListL{2,mI} = [sListL{2,mI}(1:Lsess(4)),nan,sListL{2,mI}(Lsess(4)+1:end)];
                else
                end
            end
        end
        
        
        % concentration dataset OK, 2018-10-23 AndyP
        
        nT = length(siD);
        for iT=1:nT
            mouseS(mI,iC) = mI;
            sessS(mI,iC) = iC;
            Tr(mI,iC) = iT;
            switch piD{iT}
                case {'no spot',' no spot','no spot ',' no spot '}
                    c2n = -1;
                case {'unscented', ' unscented', 'unscented ',' unscented '}
                    c2n = 0;
                case {'0.01%', '.01%',' 0.01%','0.01% ', ' 0.01% '}
                    c2n = 1;
                case {'0.1%',' 0.1%', '.1%', '0.1% ', ' 0.1% '}
                    c2n = 2;
                case {'1%',' 1%','1% ',' 1% '}
                    c2n = 3;
                case {'2%',' 2%','2% ',' 2% '}
                    c2n = 4;
                case {'nan',' nan','nan ',' nan '}
                    c2n = nan;
                otherwise
                    error('unknown concentration %d %s',iD,dateStr);
            end
            
            C(mI,iC) = c2n; % unscented->0,0.01%->1,0.1%->2,1%->3,2%->4,no spot->-1
            switch siD{iT}
                case 'A'
                    s2n = 1;
                case 'B'
                    s2n = 2;
                case 'C'
                    s2n = 3;
                case 'D'
                    s2n = 4;
                case 'E'
                    s2n = 5;
                case 'F'
                    s2n = 6;
                case 'G'
                    s2n = 7;
                case 'H'
                    s2n = 8;
                case 'I'
                    s2n = 9;
                case 'J'
                    s2n = 10;
                case 'K'
                    s2n = 11;
                case 'L'
                    s2n = 12;
                case 'M'
                    s2n = 13;
                case 'N'
                    s2n = 14;
                case {'nan',' nan','nan ',' nan '}
                    s2n = nan;
            end
            sP(mI,iC) = s2n; % sP = []; % spot position (A-N -> 1-14)
            iC = iC+1;
        end
    end
end

cx0 = nan(nM,nS);
cy0 = nan(nM,nS);
for iM=1:nM
    cd(mDir{iM});
    % get number of sessions
    %sD = dir('*.mat');
    sDprimary = sListL{1,iM};
    sDsecondary = sList{1,iM};
    
    iC = 1;
    for iS=1:nS
        
        if ~isnan(C(iM,iS))
            
            if ~isnan(sDprimary{iC})
                primaryflag = 1;
                load(sDprimary{iC},'bx','by','nx','ny','tbx','tby','ttx','tty','confb','confn','conftb','conftt');
                x0 = bx;
                y0 = by;
                nx0 = nx; 
                ny0 = ny; %#ok<NODEF>
                tx0 = ttx;
                ty0 = tty;
            else
                primaryflag = 0;
                load(sDsecondary{iC},'x0','y0','nx0','ny0','tx0','ty0');
            end
            
            % load the rest
            load(sDsecondary{iC},'cx','cy','time0','frame0','mouseLength0','trialNum','trialStartT','foundSpotT','maxSessT','radiusAroundSpot','pixpercm');
            nF = min(length(frame0),length(x0));
            if primaryflag
                x0 = x0(1:nF);
                y0 = y0(1:nF);
                nx0 = nx0(1:nF);
                ny0 = ny0(1:nF);
                tx0 = tx0(1:nF);
                ty0 = ty0(1:nF);
                tbx = tbx(1:nF);
                tby = tby(1:nF);
                confb = confb(1:nF);
                confn = confn(1:nF);
                conftb = conftb(1:nF);
                conftt = conftt(1:nF);
            else
                nx = nan(size(time0));
                tbx = nan(nF,1);
                tby = nan(nF,1);
                confb = nan(nF,1);
                confn = nan(nF,1);
                conftb = nan(nF,1);
                conftt = nan(nF,1);
            end
            
            if length(time0)>length(nx)
                %                 keyboard;
                disp(sprintf('check for misaligned video: %s',sDsecondary{iC}));
            end
            
            cx0(iM,iS) = cx;
            cy0(iM,iS) = cy;
            %             % debugging
            %             if cy < 470 && sP(iM,iS)==2
            %                 keyboard;
            %             end
            %
            w = warning('query','last');
            if ~isempty(w)
                id = w.identifier;
                warning('off',id);
            end
            if length(x0)==length(nx0)
                nP = length(x0);
            else
                nP = length(nx0);
                disp(sprintf('last nx0 points missing, truncating session %s \n',sDprimary{iS}));
            end
            x = cat(1,x,x0(1:nP));
            y = cat(1,y,y0(1:nP));
            if ~primaryflag
                usedLEAP = cat(1,usedLEAP,zeros(nP,1));
                [tnx,tny,ttx0,tty0] = switch_nose_tail(x0(1:nP),y0(1:nP),time0(1:nP),nx0,ny0,tx0(1:nP),ty0(1:nP));
                nx1 = cat(1,nx1,tnx);
                ny1 = cat(1,ny1,tny);
                tx = cat(1,tx,ttx0);
                ty = cat(1,ty,tty0);
            else
                usedLEAP = cat(1,usedLEAP,ones(nP,1));
                tnx = nx0;
                tny = ny0;
                ttx0 = tx0;
                tty0 = ty0;
                nx1 = cat(1,nx1,nx0);
                ny1 = cat(1,ny1,ny0);
                tx = cat(1,tx,tx0);
                ty = cat(1,ty,ty0);
            end
            tbx1 = cat(1,tbx1,tbx);
            tby1 = cat(1,tby1,tby);
            confb1 = cat(1,confb1,confb);
            confn1 = cat(1,confn1,confn);
            conftb1 = cat(1,conftb1,conftb);
            conftt1 = cat(1,conftt1,conftt);
            t = cat(1,t,time0(1:nP));
            sess = cat(1,sess,repmat(iS,[nP,1]));
            mouse = cat(1,mouse,repmat(iM,[nP,1]));
            trial = cat(1,trial,repmat(Tr(iM,iS),[nP,1]));
            mL = cat(1,mL,sqrt((tnx-ttx0).^2+(tny-tty0).^2)./pixpercm);
            cx1 = cat(1,cx1,repmat(cx,[nP,1]));
            cy1 = cat(1,cy1,repmat(cy,[nP,1]));
            dT = cat(1,dT,sqrt((x0(1:nP)-cx).^2+(y0(1:nP)-cy).^2)./pixpercm);
            dnT0 = sqrt((tnx-cx).^2+(tny-cy).^2)./pixpercm;
            dnT = cat(1,dnT,dnT0);
            BSx = cx-x0(1:nP);
            BSy = cy-y0(1:nP);
            BNx = tnx-x0(1:nP);
            BNy = tny-y0(1:nP);
            
            A = sqrt(BSx.^2+BSy.^2);
            B = sqrt(BNx.^2+BNy.^2);
            knonan = ~isnan(BSx) & ~isnan(BNx);
            orient0 = nan(sum(knonan),1);
            orient0(knonan) = acosd(dot([BSx(knonan),BSy(knonan)]',[BNx(knonan),BNy(knonan)]')'./(A(knonan).*B(knonan)));
            
            orient = cat(1,orient,orient0);
            dx = foaw_diff_varTs(x0(1:nP), time0(1:nP), 50, 0.2,0.1)./pixpercm;
            dy = foaw_diff_varTs(y0(1:nP), time0(1:nP), 50, 0.2,0.1)./pixpercm;
            dnx = foaw_diff_varTs(tnx, time0(1:nP), 50, 0.2,0.1)./pixpercm;
            dny = foaw_diff_varTs(tny, time0(1:nP), 50, 0.2,0.1)./pixpercm;
            bV = cat(1,bV,sqrt(dx.^2+dy.^2));
            nV = cat(1,nV,sqrt(dnx.^2+dny.^2));
            bdphi0 = zIdPhi2(dx,dy,time0(1:nP), 50, 0.2,0.1);
            ndphi0 = zIdPhi2(dnx,dny,time0(1:nP), 50, 0.2,0.1);
            bdphi = cat(1,bdphi,bdphi0);
            ndphi = cat(1,ndphi,ndphi0);
            nT = nP;
%             Dline0 = nan(nT,1);
%             theta0 = nan(nT,1);
%             slope0 = nan(nT,1);
%                         for iT=4:nT
%                             dy0 = (y0(iT)+nanmean(dy(iT-3:iT).*pixpercm))-y0(iT);
%                             dx0 = (x0(iT)+nanmean(dx(iT-3:iT).*pixpercm))-x0(iT);
%                             m0 = dy0./dx0;
%                             slope0(iT)=m0;
%                             minX = 1;
%                             maxX = 576;
%                             result = computeline([x0(iT), y0(iT)],m0, [minX maxX]);
%                             nP = length(result);
%                             x1 = result{1}(:,1);
%                             y1 = result{1}(:,2);
%                             x2 = result{nP}(:,1);
%                             y2 = result{nP}(:,2);
%             
%                             Dline0(iT) = ((y2-y1)*nx0(iT)-(x2-x1)*ny0(iT)+x2*y1-y2*x1)./sqrt((y2-y1).^2+(x2-x1).^2);
%                             dxline = sqrt((nx0(iT)-x0(iT)).^2+(ny0(iT)-x0(iT)).^2);
%                             theta0(iT) = acosd(dxline/sqrt(dxline.^2+Dline0(iT).^2));
%                         end
%             slope = cat(1,slope,slope0);
%             Dline = cat(1,Dline,Dline0./pixpercm);
%             theta = cat(1,theta,theta0);
            
            if length(foundSpotT)>length(trialStartT)
                trialNum = cat(2,trialNum,nan);
                trialStartT = cat(1,trialStartT,nan);
                flag = 1;
            else
                flag = 0;
            end
            
            nTr = cat(1,nTr,trialNum');
            f2s = cat(1,f2s,interp1(time0,frame0,foundSpotT,'nearest')-1);
            f2f = cat(1,f2f,interp1(time0,frame0,trialStartT,'nearst')-1);
            tempt = interp1(time0,frame0,foundSpotT,'nearest')-1;
            t2s = cat(1,t2s,time0(tempt));
            tempt1 = interp1(time0,frame0,trialStartT(~isnan(trialStartT)),'nearest')-1;
            if flag
                t2f = cat(1,t2f,cat(1,time0(tempt1),nan));
            else
                t2f = cat(1,t2f,time0(tempt1));
            end
            mouseT = cat(1,mouseT,repmat(iM,[length(trialNum),1]));
            sessT = cat(1,sessT,repmat(iS,[length(trialNum),1]));
            %             disp([length(trialNum),length(foundSpotT),length(trialStartT)]); %debugging
            %
            sT(iM,iS) = length(trialNum);
            if ~isempty(tempt)
                t2s1(iM,iS) = time0(tempt(1));
            else
                t2s1(iM,iS)=nan;
            end
            tdiffs = nan(length(tempt),1);
            for iPass = 1:length(tempt)
                if iPass==1
                    tdiffs(iPass)=time0(tempt(1));
                else
                    tdiffs(iPass)=time0(tempt(iPass))-time0(tempt1(iPass-1));
                end
            end
            mTlen(iM,iS) = nanmean(tdiffs);
            mTlen2(iM,iS) = nanmedian(tdiffs);
            initD(iM,iS) = dnT0(find(~isnan(dnT0),1,'first'));
            initO(iM,iS) = orient0(find(~isnan(orient0),1,'first'));
            % if length(t2s)~=length(t2f); keyboard; end % debugging
            iC = iC + 1;
            fprintf('%d/%d \n',iS,nS);
        else % do not update iC
        end
        
    end
    disp(iM);
    cd(homedir);
end

nx = nx1;
ny = ny1;
tbx = tbx1;
tby = tby1;
confb = confb1;
confn = confn1;
conftb = conftb1;
conftt = conftt1;
edge = ~(x>10*5.1744 & y>10*5.1744 & x<576-10*5.1744 & y<480-10*5.1744 & sqrt((x-576).^2+(y-1).^2)>100); % 6cm from edge (@5.1744 pixels / cm), exclude feeder zone in circle
end