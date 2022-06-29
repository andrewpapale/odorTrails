function [trainFeatures, testFeatures,featureindices] = helperExtractFeatures(trainData,testData,T,AR_order,level)
% This function is only in support of XpwWaveletMLExample. It may change or
% be removed in a future release.
trainFeatures = [];
testFeatures = [];

for idx =1:size(trainData,1)
    x = trainData(idx,:);
    k = find(~isnan(x));
    x(k) = detrend(x(k),0,'Continuous',false);
    arcoefs = blockAR(x,AR_order,T);
    se = shannonEntropy(x,T,level);
    [cp,rh] = leaders(x,T);
    wvar = modwtvar(modwt(x(k),'db2'),'db2');
    trainFeatures = [trainFeatures; arcoefs se cp rh wvar']; %#ok<AGROW>
    
end

for idx =1:size(testData,1)
    x1 = testData(idx,:);
    k = find(~isnan(x1));
    x1(k) = detrend(x1(k),0);
    arcoefs = blockAR(x1,AR_order,T);
    se = shannonEntropy(x1,T,level);
    [cp,rh] = leaders(x1,T);
    wvar = modwtvar(modwt(x1(k),'db2'),'db2');
    testFeatures = [testFeatures;arcoefs se cp rh wvar']; %#ok<AGROW>
    
end

featureindices = struct();
% 4*8
featureindices.ARfeatures = 1:32;
startidx = 33;
endidx = 33+(16*8)-1;
featureindices.SEfeatures = startidx:endidx;
startidx = endidx+1;
endidx = startidx+7;
featureindices.CP2features = startidx:endidx;
startidx = endidx+1;
endidx = startidx+7;
featureindices.HRfeatures = startidx:endidx;
startidx = endidx+1;
endidx = startidx+13;
featureindices.WVARfeatures = startidx:endidx;
end


function se = shannonEntropy(x,numbuffer,level)
numwindows = numel(x)/numbuffer;
y = buffer(x,numbuffer);
[i,j]=find(~isnan(y));
se = zeros(2^level,size(y,2));
for kk = 1:size(y,2)
    if ~isempty(i(j==kk))
        wpt = modwpt(y(i(j==kk),kk),level);
        % Sum across time
        E = sum(wpt.^2,2);
        Pij = wpt.^2./E;
        % The following is eps(1)
        se(:,kk) = -sum(Pij.*log(Pij+eps),2);
    end
end
se = reshape(se,2^level*numwindows,1);
se = se';
end


function arcfs = blockAR(x,order,numbuffer)
numwindows = numel(x)/numbuffer;
y = buffer(x,numbuffer);
[i,j] = find(~isnan(y));
arcfs = zeros(order,size(y,2));
for kk = 1:size(y,2)
    if ~isempty(i(j==kk))
        artmp =  arburg(y(i(j==kk),kk),order);
        arcfs(:,kk) = artmp(2:end);
    end
end
arcfs = reshape(arcfs,order*numwindows,1);
arcfs = arcfs';
end


function [cp,rh] = leaders(x,numbuffer)
y = buffer(x,numbuffer);
[i,j]=find(~isnan(y));
cp = zeros(1,size(y,2));
rh = zeros(1,size(y,2));
for kk = 1:size(y,2)
    if ~isempty(i(j==kk))
        [~,h,cptmp] = dwtleader(y(i(j==kk),kk));
        cp(kk) = cptmp(2);
        rh(kk) = range(h);
    end
end
end


