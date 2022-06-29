% Wrap_RandWalk

nM = max(mouse);
nS = max(sess);

x_sim = [];
y_sim = [];
th = [];
dnT_sim = [];
idphi_sim = [];
nV_sim = [];
nV1_sim = [];
frame_sim = [];
mouse_sim = [];
sess_sim = [];
for iM=1:nM
    disp(iM);
    for iS=1:nS
        k = mouse==iM & sess==iS;
        kT = mouseT==iM & sessT==iS;
        
        if sum(k)>0
            
            xT=nanmedian(xT1);
            yT=nanmedian(yT1);
            
            frame0 = frame(k);
            vdist = nV(k);
            
            [x0,y0,th0,dnT0,idphi0,v0,V0]=randwalk_01(frame0,vdist,[xT,yT]);
            
            frame_sim = cat(1,frame_sim,(1:length(x0))');
            mouse_sim = cat(1,mouse_sim,repmat(iM,[length(x0),1]));
            sess_sim = cat(1,sess_sim,repmat(iS,[length(x0),1]));
            x_sim = cat(1,x_sim,x0);
            y_sim = cat(1,y_sim,y0);
            th = cat(1,th,th0);
            dnT_sim = cat(1,dnT_sim,dnT0);
            idphi_sim = cat(1,idphi_sim,idphi0);
            nV_sim = cat(1,nV_sim,v0);
            nV1_sim = cat(1,nV1_sim,V0);
            
        end
    end
end