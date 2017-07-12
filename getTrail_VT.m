function [xT,yT] = getTrail_VT(arena_data)
% 2017-06-09 AndyP
% [xT,yT] = getTrail_VT(arena_data);
% get coordinates of trail (xT, yT) based on interpolation from manual
% clicking
% when finished, press the 'shift' key

mI = double(arena_data.MedianImage);

figure(1); clf; imagesc(mI); 

% 
% iT=1;
% flag = 0;
% while ~flag
%     [x0,y0] = ginput(1);
%     if ~isempty(x0)
%         xT(iT)=x0;
%         yT(iT)=y0;
%         key = get(gcf,'CurrentKey');
%         plot(xT(iT),yT(iT),'r.','markersize',10); hold on;
%         iT=iT+1;
%     end
%     if(strcmp (key , 'shift'))
%         flag = 1;
%     end
% end
% 
% nameStr = strsplit(arena_data.videofilename,'\');
% endSplit = strsplit(nameStr,'.');
% fileName = strcat(endSplit{1},'_Trail_Coordinates');
% 
% save(fileName,'xT','yT');

%end



% doFill = false;
% tolerance = 100;

% remove edges
mI(1:50,:)=nan;
mI(:,1:50)=nan;
mI(end-50:end,:)=nan;
mI(:,end-50:end)=nan;

figure(2); clf;  imagesc(trail0);


BW2 = edge(mI,'canny',[minThr,maxThr],sigma);

trail0 = medfilt2(trail0);
trail0 = bwareaopen(trail0,40);

% if doFill  % connect closest two points throughout image
%     [y,x]=find(trail0==1);
%     bins = linspace(0,100,100);
%     iC=1;
%     x0 = []; y0 = []; ix = [];
%     for i=1:length(y)
%         % find distance btw x(i), y(i) and all other points...
%         D0 = sqrt((x-x(i)).^2+(y-y(i)).^2);
%         H = hist(D0,bins);
%         if H(2)<4 % find edge points ~ those that border < 4 points..
%             D(iC,:)=D0;
%             x0 = cat(1,x0,x(i));
%             y0 = cat(1,y0,y(i));
%             ix = cat(1,ix,i);
%             iC=iC+1;
%         end
%     end
%     nC = length(ix); % number of edge points
%     % get subset of distances between edge points
%     for iC=1:nC
%         for iD=1:nC
%             D1(iC,iD)=D(iC,ix(iD));
%         end
%     end
%     % connect 2 edge points closer together than tolerance
%     D1(D1==0)=tolerance+1;
%     imagesc(trail0); hold on;
%     iF = 1;
%     for iC=1:nC
%         flag = 0;
%         [D2,newI] = sort(D1(iC,:)); % sort in order of increasing distance
%         for iD=1:nC
%             if D2(1,iD)<tolerance && D2(1,iD)>1 && ~flag
%                 %flag = 1;  % only connect one point per edge
% %                 Pc = plot(x0(iC),y0(iC),'r.','markersize',20); hold on;
% %                 Pd = plot(x0(newI(iD)),y0(newI(iD)),'g.','markersize',20); hold on;
%                 a = [x0(iC),y0(iC)];
%                 b = [x0(newI(iD)),y0(newI(iD))];
%                 %keyboard;
%                 result = linepts(a, b);
%                 for iT=1:length(result);
%                         fill{iF}(iT,:)=result{iT}; 
%                 end 
%                 iF = iF + 1;
%             end
%         end
%     end
    
    figure(3); clf;  imagesc(trail0);
    
%     nF = length(fill);
%     for iF=1:nF
%         fillx = fill{iF}(:,1);
%         filly = fill{iF}(:,2);
%         nX0 = length(fillx);
%         if ~any(isnan(fillx(:))) && ~any(isnan(filly(:)))
%             for iX=1:nX0
%                 trail0(filly(iX),fillx(iX))=1;
%             end
%         end
%     end
%     
    
end

figure(4); clf;  imagesc(trail0);

[yT,xT] = find(trail0==1);

% keyboard;
% iT=1;
% flag = 0;
% while ~flag
%     [x0,y0] = ginput(1);
%     if ~isempty(x0)
%         xT0(iT)=x0;
%         yT0(iT)=y0;
%         key = get(gcf,'CurrentKey');
%         
%         
%         
%         
%         plot(xT(iT),yT(iT),'r.','markersize',10); hold on;
%         iT=iT+1;
%     end
%     if(strcmp (key , 'shift'))
%         flag = 1;
%     end
% end



end


function result = linepts(a, b)

m = (a(2) - b(2)) / (a(1) - b(1));
n = b(2) - b(1) * m;

x = min(a(1), b(1)) : max(a(1), b(1));
y = m * x + n;
y = round(y);

for i = 1 : length(x)
    result{i} = [x(i), y(i)];
end

end

