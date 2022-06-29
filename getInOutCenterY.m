
doAll = true;
Dpix = 112;
tol = 10;
timeW = 50;

nM = max(mouse);
nS = max(sess);
angIn = [];
angOut = [];
thetaA = [];
thetaB = [];
thetaC = [];
ABCIn = [];
minDIn = [];
ABCOut = [];
minDOut = [];
dTA = [];
dTB = [];
dTC = [];
mouse0 = [];
sess0 = [];
count0 = [];
iC = 1;
for iM=1:nM
    disp(iM);
    for iS=1:nS
        disp(iS);
        k = mouse==iM & sess==iS;
        kT = mouseT==iM & sessT==iS;
        kR = dTc < Dpix/11.2;
        
        if sum(k)>0
            %F = figure(1); clf;
            %plot(xT1(kT),yT1(kT),'r.');
           % hold on;
            
            %         scatter(x(k),y(k),5,dTc(k));
            %         caxis([0 100/11.2]);
            %        
            
           % plot(x(k),y(k),'k.','markersize',1);
            
            
            %        pause;
            
            % get entry angle
            
            tR = round(smooth(+(k & kR),timeW,'moving'));
            tin = cat(1,0,diff(tR))>0;
            tout = cat(1,0,diff(tR))<0;
            tin = find(tin==1);
            tout = find(tout==1);
            
            k0 = tout-tin >= timeW;
            
            
            
            %plot(x(tin),y(tin),'bx','markersize',20);
            %plot(x(tout),y(tout),'gx','markersize',20);
            %disp(sum(tin));
            %disp(sum(tout));
            %pause;
            % to do...check that ABCIn and ABCOut correspond to closest
            % dTarm
            
            dTA0 = getNearestNonNan(dTarm{1},tin,tol);
            dTB0 = getNearestNonNan(dTarm{2},tin,tol);
            dTC0 = getNearestNonNan(dTarm{3},tin,tol);
            kA = dTA0-dTB0 > 0 & dTA0-dTC0 > 0;
            kB = dTB0-dTA0 > 0 & dTB0-dTC0 > 0;
            kC = dTC0-dTA0 > 0 & dTC0-dTB0 > 0;
            
            dTA0 = getNearestNonNan(dTarm{1},tout,tol);
            dTB0 = getNearestNonNan(dTarm{2},tout,tol);
            dTC0 = getNearestNonNan(dTarm{3},tout,tol);
            kE = dTA0-dTB0 > 0 & dTA0-dTC0 > 0;
            kF = dTB0-dTA0 > 0 & dTB0-dTC0 > 0;
            kG = dTC0-dTA0 > 0 & dTC0-dTB0 > 0;
            
            k0 = ~isnan(dTA0) & k0;
            if ~doAll
                k0 = (kA | kB | kC) & (kE | kF | kG) & k0;
            elseif doAll
                
            end
            %disp(sum(k0));
            
%             for iT=1:length(tin);
%                 if k0(iT)
%                     x0 = x(k);
%                     y0 = y(k);
%                     T2 = plot(x0(tout(iT):end),y0(tout(iT):end),'b.','markersize',2);
%                     T1 = scatter(nearX(tin(iT):tout(iT)),nearY(tin(iT):tout(iT)),5,frame(tin(iT):tout(iT)),'filled');
%                     %                 plot(nearX(tin(iT)),nearY(tin(iT)),'b.','markersize',20);
%                     %                 plot(nearX(tout(iT)),nearY(tout(iT)),'b.','markersize',20);
%                     pause;
%                     delete(T1);
%                     delete(T2);
%                 end
%             end
            
            dTA0 = getNearestNonNan(dTarm{1},tin,tol);
            dTB0 = getNearestNonNan(dTarm{2},tin,tol);
            dTC0 = getNearestNonNan(dTarm{3},tin,tol);
            [minD0,ABCIn0]=nanmin(cat(2,dTA0(k0==1),dTB0(k0==1),dTC0(k0==1)),[],2);
            minDIn = cat(1,minDIn,minD0);
            ABCIn = cat(1,ABCIn,ABCIn0);
            
            %             disp(ABCIn0);
            %             pause;
            
            dTA0 = getNearestNonNan(dTarm{1},tout,tol);
            dTB0 = getNearestNonNan(dTarm{2},tout,tol);
            dTC0 = getNearestNonNan(dTarm{3},tout,tol);
            [minD0,ABCOut0]=nanmin(cat(2,dTA0(k0==1),dTB0(k0==1),dTC0(k0==1)),[],2);
            minDOut = cat(1,minDOut,minD0);
            ABCOut = cat(1,ABCOut,ABCOut0);
            
            %             for iL=1:length(k0)
            %                 if ABCIn(iL)==1 & ABCOut(iL)==1
            %                     keyboard;
            %                 end
            %             end
            
            angA = acos(getNearestNonNan(dotA,tin,tol))*180/pi;
            angB = acos(getNearestNonNan(dotB,tin,tol))*180/pi;
            angC = acos(getNearestNonNan(dotC,tin,tol))*180/pi;
            AngAll = cat(2,angA(k0==1),angB(k0==1),angC(k0==1));
            angIn = cat(1,angIn,AngAll(ABCIn0));
            
            angA = acos(getNearestNonNan(dotA,tout,tol))*180/pi;
            angB = acos(getNearestNonNan(dotB,tout,tol))*180/pi;
            angC = acos(getNearestNonNan(dotC,tout,tol))*180/pi;
            AngAll = cat(2,angA(k0==1),angB(k0==1),angC(k0==1));
            angOut = cat(1,angOut,AngAll(ABCOut0));
            
            thetaA = cat(1,thetaA,repmat(mode(theta1(k)),[length(ABCOut0),1]));
            thetaB = cat(1,thetaB,repmat(mode(theta2(k)),[length(ABCOut0),1]));
            thetaC = cat(1,thetaC,repmat(mode(theta3(k)),[length(ABCOut0),1]));
            
            All = (mode(theta1(k))+mode(theta2(k))+mode(theta3(k)));
            
%             if All<300
%                 keyboard;   % bent inward y #51
%             end
            
            mouse0 = cat(1,mouse0,repmat(iM,[length(ABCOut0),1]));
            sess0 = cat(1,sess0,repmat(iS,[length(ABCOut0),1]));
            count0 = cat(1,count0,repmat(iC,[length(ABCOut0),1]));
            
            iC=iC+1;
            
        end
    end
end
