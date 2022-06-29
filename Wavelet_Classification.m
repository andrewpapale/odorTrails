% find max time per lap
maxt = 1;
for iM=1:nM
    for iS=1:nS
        k = mouse==iM & sess==iS;
        kT = mouseT==iM & sessT==iS;
        L = nTr(kT);
        nT = sum(kT);
        for iL=1:nT
            if sum(t2s0(k)==iL)>maxt
                maxt = nansum(t2s0(k));
            end
        end
    end
end

% get laps in nL x nSamples format.  Buffer ends with nans
nTL = size(mouseT,1);
Data = nan(nTL,maxt);
iC = 1;
for iM=1:nM
    for iS=1:nS
        k = mouse==iM & sess==iS;
        kT = mouseT==iM & sessT==iS;
        L = nTr(kT);
        nT = sum(kT);
        for iT=1:min(nT,5)
            L0 = L(iT);
            k0 = k & t2s0==iT & confb > 0.5 & confn > 0.5;
            Data(iC,1:sum(k0))=x(k0);
            %LabelConc{iC} = mat2str(C(iM,iS));
            %LabelLap{iC} = mat2str(iT);
            LabelSp{iC} = mat2str(sP(iM,iS));
            iC = iC + 1;
        end
    end
    disp(iM);
end

timeWindow = floor(maxt./10);
AFData.Data = Data(:,1:10*timeWindow);
%AFData.Labels = LabelConc';
AFData.Labels = LabelSp';

percent_train = 70;
ARorder = 4;
MODWPTlevel = 5;


[trainData,testData,trainLabels,testLabels,trainIdx,testIdx] = helperRandomSplit1(percent_train,AFData);
[trainFeatures,testFeatures,featureindices] = ...
    helperExtractFeatures1(trainData,testData,timeWindow,ARorder,MODWPTlevel);
features = [trainFeatures; testFeatures];


nT = size(features,1);
skip = zeros(nT,1);
train = zeros(nT,1);
test = zeros(nT,1);
for iT=1:nT
    if any(trainIdx==iT)
        train(iT) =1;
    else
        test(iT)=1;
    end
    if ~all(isnan(features(iT,:)))
        skip(iT) = 0;
    else
        skip(iT) = 1;
    end
end

iC = 1; iD = 1;
for iT=1:nT
    if train(iT) && skip(iT)
        trainLabels(iC)={nan};
        iC = iC+1;
    end
    if test(iT) && skip(iT)
        testLabels(iD)={nan};
        iD = iD+1;
    end
end

trainLabels(cellfun(@isnan,trainLabels))=[];
testLabels(cellfun(@isnan,testLabels))=[];
features(skip==1,:)=[];


rng(1)
template = templateSVM(...
    'KernelFunction','polynomial',...
    'PolynomialOrder',2,...
    'KernelScale','auto',...
    'BoxConstraint',1,...
    'Standardize',true);
model = fitcecoc(...
    features,...
    [trainLabels;testLabels],...
    'Learners',template,...
    'Coding','onevsone',...
    'ClassNames',{'1','2','3','4','5','6','7','8','9','10','11','12','13','14'});
kfoldmodel = crossval(model,'KFold',5);
classLabels = kfoldPredict(kfoldmodel);
loss = kfoldLoss(kfoldmodel)*100
[confmatCV,grouporder] = confusionmat([trainLabels;testLabels],classLabels);

