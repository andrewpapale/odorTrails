function [GeoDist]=GeodesicDistanceBetweenTwoPointsOnImage(Im,x0,y0,x1,y1,WeightDist,NumSteps)
%Find minimal Geodesic distance between point two point (x0,y0 and x1,y1) on image (Im)
%The image can be color or grayscale
%return geodesic distance between x0,y0 and x1,y1

%Use the image as topological map with color as height and use the
%Dijkstra's algorithm to find minimal geodesic distance between x1,y1 and x0,y0

%Input:

%Im color image (or other multichannel image)

%x0,y0 coordinates (on the image) of the origin point
%x1,y1 Cordinates of the target point

%NumSteps optional parameter that allow  you to limit the number of cycle of the calculation this may result faster but less accurate calculation

%WeightDist the distance component have two parts distance travel on
%the image plane and distance in the color/grayscale value this parameter
%control the relative weight of the distance on the image plane

%Output:
%DistMap a map of geodesic distance to coordinate x0,y0

%Method:
%Use  the image as topological map and use the Dijkstra's algorithm to find geodesic distance


if nargin < 4
    WeightDist=1;
end

if nargin < 5
    NumSteps=1000;
end


[Sy,Sx,~]=size(Im);% size of image
DistMap=ones(Sy,Sx)*10000000;%Geodesic distance map from x0,y0
DistMap(y0,x0)=0;%Distance to origin is zero
PrevCycleP=zeros(10000,2);% Coordinates of points their Geodesic distance was updated in previous cycle
NPrev=1; %Number of point their geodesic distance was updated in previous cycle
ThisCycleP=zeros(10000,2);%
PrevCycleP(1,:)=[x0,y0];%The initial point for the scan
for f=1:NumSteps% Propogate distance x0,y0 in cycles using the Dijkstra's algorithm
    %%----------------------Create geodesic  distance map  ----------------------------------------------------
    NThis=0;
    for fPr=1:NPrev % Find the distance from each point wh that have its distance calculated in previous
        x=PrevCycleP(fPr,1);% Coordinate of point from previous cycle
        y=PrevCycleP(fPr,2);% Coordinate of point from previous cycle
        if (~(x>0 && y>0 && x<=Sx && y<=Sy)) % check index is in image limits
            continue;
        end
        
        for nx=x-1:x+1%Scan all direct neighbours of  current point (x,y)
            for ny=y-1:y+1%Scan all direct neighbours of  current point
                if (nx>0 && ny>0 && nx<=Sx && ny<=Sy && ~(ny==y && nx==x)) % See that the neighbour is within image limits
                    Dst=sqrt(sum((Im(y,x,:)-Im(ny,nx,:)).^2))+WeightDist*sqrt((nx-x)^2+(ny-y)^2)+DistMap(y,x);% Distance to neighbours points color distace to neighbour+ distance on plane+ distance to origin point
                    if (Dst<DistMap(ny,nx))% Update point if the distance shorter then previously found distance
                        DistMap(ny,nx)=Dst;% Update new distance for the point
                        NThis=NThis+1;% Number of points their distance was updated this cycle
                        ThisCycleP(NThis,:)=[nx,ny];% Add point to the list of point whose distance was updated this cycle
                    end
                end
            end
        end
    end
    %..............Update parameters for next cycle..........................................................................................
    NPrev=NThis;%Update number of starting points for next cycle
    [PrevCycleP, ThisCycleP] = deal(ThisCycleP, PrevCycleP);%Update coordinate of startin point for next cycle
    %display(f);
    if NThis==0
        break;
    end
end
GeoDist=DistMap(y1,x1);
%figure, imshow(DistMap/max(max(DistMap)));% Display distance map
end

            
       