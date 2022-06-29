function [trainData, testData, trainLabels, testLabels] = helperRandomSplit(percent_train_split,ECGData)
% This function is only in support of XpwWaveletMLExample. It may change or
% be removed in a future release.
    Labels = ECGData.Labels;
    Data = ECGData.Data;
    percent_train_split = percent_train_split/100;
    idxARRbegin = find(strcmpi(Labels,'ARR'),1,'first');
    idxARRend = find(strcmpi(Labels,'ARR'),1,'last');
    Narr = idxARRend-idxARRbegin+1;
    idxCHFbegin = find(strcmpi(Labels,'CHF'),1,'first');
    idxCHFend = find(strcmpi(Labels,'CHF'),1,'last');
    Nchf = idxCHFend-idxCHFbegin+1;
    idxNSRbegin = find(strcmpi(Labels,'NSR'),1,'first');
    idxNSRend = find(strcmpi(Labels,'NSR'),1,'last');
    Nnsr = idxNSRend-idxNSRbegin+1;
    % Obtain number needed for percentage split
    num_train_arr = round(percent_train_split*Narr);
    num_train_chf = round(percent_train_split*Nchf);
    num_train_nsr = round(percent_train_split*Nnsr);
    rng default;
    Parr = randperm(Narr,num_train_arr);
    Pchf = randperm(Nchf,num_train_chf);
    Pnsr = randperm(Nnsr,num_train_nsr);
    notParr = setdiff(1:Narr,Parr);
    notPchf = setdiff(1:Nchf,Pchf);
    notPnsr = setdiff(1:Nnsr,Pnsr);
    ARRdata = Data(idxARRbegin:idxARRend,:);
    ARRLabels = Labels(idxARRbegin:idxARRend);
    CHFdata = Data(idxCHFbegin:idxCHFend,:);
    CHFLabels = Labels(idxCHFbegin:idxCHFend);
    NSRdata = Data(idxNSRbegin:idxNSRend,:);
    NSRLabels = Labels(idxNSRbegin:idxNSRend);
    trainARR = ARRdata(Parr,:);
    trainARRLabels = ARRLabels(Parr);
    testARR = ARRdata(notParr,:);
    testARRLabels = ARRLabels(notParr);
    trainCHF = CHFdata(Pchf,:);
    trainCHFLabels = CHFLabels(Pchf);
    testCHF = CHFdata(notPchf,:);
    testCHFLabels = CHFLabels(notPchf);
    trainNSR = NSRdata(Pnsr,:);
    trainNSRLabels = NSRLabels(Pnsr);
    testNSR = NSRdata(notPnsr,:);
    testNSRLabels = NSRLabels(notPnsr);
    trainData = [trainARR ; trainCHF; trainNSR];
    trainLabels = [trainARRLabels ; trainCHFLabels; trainNSRLabels];
    testData = [testARR ; testCHF; testNSR];
    testLabels = [testARRLabels; testCHFLabels; testNSRLabels];
   
    