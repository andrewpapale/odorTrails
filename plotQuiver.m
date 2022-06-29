subplot(4,2,1); 
quiver(X,Y,cos(Htop(:,:,2)'),sin(Htop(:,:,2))'); 
viscircles([cx0(wP==8 & ~removeS),cy0(wP==8 & ~removeS)],repmat(5.174*1.5,[sum(wP(~removeS)==8),1]));

subplot(4,2,2); 
quiver(X,Y,cos(Htop(:,:,1)'),sin(Htop(:,:,1))'); 
viscircles([cx0(wP==1 & ~removeS),cy0(wP==1 & ~removeS)],repmat(5.174*1.5,[sum(wP(~removeS)==1),1]));

subplot(4,2,3); 
quiver(X,Y,cos(Hrig(:,:,2)'),sin(Hrig(:,:,2))'); 
viscircles([cx0(wP==2 & ~removeS),cy0(wP==2 & ~removeS)],repmat(5.174*1.5,[sum(wP(~removeS)==2),1]));

subplot(4,2,4); 
quiver(X,Y,cos(Hrig(:,:,1)'),sin(Hrig(:,:,1))'); 
viscircles([cx0(wP==3 & ~removeS),cy0(wP==3 & ~removeS)],repmat(5.174*1.5,[sum(wP(~removeS)==3),1]));

subplot(4,2,5); 
quiver(X,Y,cos(Hbot(:,:,2)'),sin(Hrig(:,:,2))'); 
viscircles([cx0(wP==5 & ~removeS),cy0(wP==5 & ~removeS)],repmat(5.174*1.5,[sum(wP(~removeS)==5),1]));

subplot(4,2,6); 
quiver(X,Y,cos(Hbot(:,:,1)'),sin(Hrig(:,:,1))'); 
viscircles([cx0(wP==4 & ~removeS),cy0(wP==4 & ~removeS)],repmat(5.174*1.5,[sum(wP(~removeS)==4),1]));

subplot(4,2,7); 
quiver(X,Y,cos(Hlef(:,:,2)'),sin(Hlef(:,:,2))'); 
viscircles([cx0(wP==7 & ~removeS),cy0(wP==7 & ~removeS)],repmat(5.174*1.5,[sum(wP(~removeS)==7),1]));

subplot(4,2,8); 
quiver(X,Y,cos(Hlef(:,:,1)'),sin(Hlef(:,:,1))'); 
viscircles([cx0(wP==6 & ~removeS),cy0(wP==6 & ~removeS)],repmat(5.174*1.5,[sum(wP(~removeS)==6),1]));

