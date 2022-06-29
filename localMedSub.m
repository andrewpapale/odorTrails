function tempF = localMedSub(temp)



[X,Y]=meshgrid((1:1920),(1:1080));
temp = adapthisteq(temp); 
temp = double(temp);
tempF = nan(size(temp));
nX = 1:100:1920;
nY = 1:100:1080;
for iX=1:length(nX)
    for iY=1:length(nY)
        if iX<length(nX) && iY < length(nY)
            k = X >= nX(iX) & X < nX(iX+1) & Y >= nY(iY) & Y < nY(iY+1);
            type = 1;
        elseif iX==length(nX) && iY < length(nY)
            k = X >= nX(iX) & Y >= nY(iY) & Y < nY(iY+1);
            type = 2;
        elseif iX<length(nX) && iY== length(nY)
            k = X >= nX(iX) & X < nX(iX+1) & Y >= nY(iY);
            type = 3;
        elseif iX==length(nX) && iY==length(nY)
            k = X>=nX(iX) & Y>=nY(iY);
            type = 4;
        end
        switch type
            case 1
                Xb = 100;
                Yb = 100;
            case 2
                Xb = 20;
                Yb = 100;
            case 3
                Xb = 100;
                Yb = 80;
            case 4
                Xb = 20;
                Yb = 80;
        end
        ii=find(k(:)==1);
        [ix,iy] = ind2sub([1080,1920],ii);
        Fij = fit([X(k),Y(k)],temp(k),'poly55');
        Z = feval(Fij,[X(k),Y(k)]);
        temp0 = temp(k)-Z;
        tempF(k)=reshape(temp0,[Xb,Yb]);
    end
end