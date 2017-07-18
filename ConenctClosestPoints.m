function data = ConenctClosestPoints(data)

%
% connect closest components
%
Skeleton =bwmorph(data,'skeleton');
endpoints = bwmorph(Skeleton,'endpoints');
[y0,x0] = find(endpoints==1);
nC = length(x0); % number of edge points
% get distances between endpoints points
D1 = nan(nC,nC);
for iC=1:nC
    for iD=1:nC
        D1(iC,iD)= sqrt((x0(iC)-x0(iD)).^2+(y0(iC)-y0(iD)).^2);
    end
end
D1(D1==0)=nan;
% connect each endpoint to its closest neighbor based on D1

for iC=1:nC
    [~,index] = nanmin(D1(iC,:));
    a = [x0(iC),y0(iC)];
    b = [x0(index),y0(index)];
    result = linepts(a, b);
    fill = [];
    for iT=1:length(result);
        fill(iT,:)=result{iT}; %#ok<AGROW>
        if ~isnan(fill(iT,2)) && ~isnan(fill(iT,1))
            data(fill(iT,2),fill(iT,1))=1;
        end
    end
%     F = figure(1); clf; imagesc(data); hold on;
%     Pc = plot(x0(iC),y0(iC),'r.','markersize',20); hold on;
%     Pd = plot(x0(index),y0(index),'g.','markersize',20); hold on;
%     plot(fill(:,1),fill(:,2),'bx','markersize',10);
end

end