function [trainData, testData, trainLabels, testLabels,trainIdx,testIdx] = helperRandomSplit1(percent_train_split,ECGData)
% This function is only in support of XpwWaveletMLExample. It may change or
% be removed in a future release.
Labels = ECGData.Labels;
Data = ECGData.Data;
percent_train_split = percent_train_split/100;
idxA = find(strcmpi(Labels,'1'));
Na = length(idxA);
idxB = find(strcmpi(Labels,'2'));
Nb = length(idxB);
idxC = find(strcmpi(Labels,'3'));
Nc = length(idxC);
idxD = find(strcmpi(Labels,'4'));
Nd = length(idxD);
idxE = find(strcmpi(Labels,'5'));
Ne = length(idxE);
idxF = find(strcmpi(Labels,'5'));
Nf = length(idxF);
idxG = find(strcmpi(Labels,'5'));
Ng = length(idxG);
idxH = find(strcmpi(Labels,'5'));
Nh = length(idxH);
idxI = find(strcmpi(Labels,'5'));
Ni = length(idxI);
idxJ = find(strcmpi(Labels,'5'));
Nj = length(idxJ);
idxK = find(strcmpi(Labels,'5'));
Nk = length(idxK);
idxL = find(strcmpi(Labels,'5'));
Nl = length(idxL);
idxM = find(strcmpi(Labels,'5'));
Nm = length(idxM);
idxN = find(strcmpi(Labels,'5'));
Nn = length(idxN);


% Obtain number needed for percentage split
num_train_a = round(percent_train_split*Na);
num_train_b = round(percent_train_split*Nb);
num_train_c = round(percent_train_split*Nc);
num_train_d = round(percent_train_split*Nd);
num_train_e = round(percent_train_split*Ne);
num_train_f = round(percent_train_split*Nf);
num_train_g = round(percent_train_split*Ng);
num_train_h = round(percent_train_split*Nh);
num_train_i = round(percent_train_split*Ni);
num_train_j = round(percent_train_split*Nj);
num_train_k = round(percent_train_split*Nk);
num_train_l = round(percent_train_split*Nl);
num_train_m = round(percent_train_split*Nm);
num_train_n = round(percent_train_split*Nn);
rng default;
Pa = randperm(Na,num_train_a);
Pb = randperm(Nb,num_train_b);
Pc = randperm(Nc,num_train_c);
Pd = randperm(Nd,num_train_d);
Pe = randperm(Ne,num_train_e);
Pf = randperm(Ne,num_train_f);
Pg = randperm(Ne,num_train_g);
Ph = randperm(Ne,num_train_h);
Pi = randperm(Ne,num_train_i);
Pj = randperm(Ne,num_train_j);
Pk = randperm(Ne,num_train_k);
Pl = randperm(Ne,num_train_l);
Pm = randperm(Ne,num_train_m);
Pn = randperm(Ne,num_train_n);
notPa = setdiff(1:Na,Pa);
notPb = setdiff(1:Nb,Pb);
notPc = setdiff(1:Nc,Pc);
notPd = setdiff(1:Nd,Pd);
notPe = setdiff(1:Ne,Pe);
notPf = setdiff(1:Nf,Pf);
notPg = setdiff(1:Ng,Pg);
notPh = setdiff(1:Nh,Ph);
notPi = setdiff(1:Ni,Pi);
notPj = setdiff(1:Nj,Pj);
notPk = setdiff(1:Nk,Pk);
notPl = setdiff(1:Nl,Pl);
notPm = setdiff(1:Nm,Pm);
notPn = setdiff(1:Nn,Pn);
Adata = Data(idxA,:);
ALabels = Labels(idxA);
Bdata = Data(idxB,:);
BLabels = Labels(idxB);
Cdata = Data(idxC,:);
CLabels = Labels(idxC);
Ddata = Data(idxD,:);
DLabels = Labels(idxD);
Edata = Data(idxE,:);
ELabels = Labels(idxE);
Fdata = Data(idxF,:);
FLabels = Labels(idxF);
Gdata = Data(idxG,:);
GLabels = Labels(idxG);
Hdata = Data(idxH,:);
HLabels = Labels(idxH);
Idata = Data(idxI,:);
ILabels = Labels(idxI);
Jdata = Data(idxJ,:);
JLabels = Labels(idxJ);
Kdata = Data(idxK,:);
KLabels = Labels(idxK);
Ldata = Data(idxL,:);
LLabels = Labels(idxL);
Mdata = Data(idxM,:);
MLabels = Labels(idxM);
Ndata = Data(idxN,:);
NLabels = Labels(idxN);
trainA = Adata(Pa,:);
trainALabels = ALabels(Pa);
testA = Adata(notPa,:);
testALabels = ALabels(notPa);
trainB = Bdata(Pb,:);
trainBLabels = BLabels(Pb);
testB = Bdata(notPb,:);
testBLabels = BLabels(notPb);
trainC = Cdata(Pc,:);
trainCLabels = CLabels(Pc);
testC = Cdata(notPc,:);
testCLabels = CLabels(notPc);
trainD = Ddata(Pd,:);
trainDLabels = DLabels(Pd);
testD = Ddata(notPd,:);
testDLabels = DLabels(notPd);
trainE = Edata(Pe,:);
trainELabels = ELabels(Pe);
testE = Edata(notPe,:);
testELabels = ELabels(notPe);
trainF = Fdata(Pf,:);
trainFLabels = FLabels(Pf);
testF = Fdata(notPf,:);
testFLabels = FLabels(notPf);
trainG = Gdata(Pg,:);
trainGLabels = GLabels(Pg);
testG = Gdata(notPg,:);
testGLabels = GLabels(notPg);
trainH = Hdata(Ph,:);
trainHLabels = HLabels(Ph);
testH = Hdata(notPh,:);
testHLabels = HLabels(notPh);
trainI = Idata(Pi,:);
trainILabels = ILabels(Pi);
testI = Idata(notPi,:);
testILabels = ILabels(notPi);
trainJ = Jdata(Pj,:);
trainJLabels = JLabels(Pj);
testJ = Jdata(notPj,:);
testJLabels = JLabels(notPj);
trainK = Kdata(Pk,:);
trainKLabels = KLabels(Pk);
testK = Kdata(notPk,:);
testKLabels = KLabels(notPk);
trainL = Ldata(Pl,:);
trainLLabels = LLabels(Pl);
testL = Ldata(notPl,:);
testLLabels = LLabels(notPl);
trainM = Mdata(Pm,:);
trainMLabels = MLabels(Pm);
testM = Mdata(notPm,:);
testMLabels = MLabels(notPm);
trainN = Ndata(Pn,:);
trainNLabels = NLabels(Pn);
testN = Ndata(notPn,:);
testNLabels = NLabels(notPn);



trainIdx = [Pa,Pb,Pc,Pd,Pe,Pf,Pg,Ph,Pi,Pj,Pk,Pl,Pm,Pn];
testIdx = [notPa,notPb,notPc,notPd,notPe,notPf,notPg,notPh,notPi,notPj,notPk,notPl,notPm,notPn];
trainData = [trainA ; trainB; trainC; trainD; trainE; trainF; trainG; trainH; trainI; trainJ; trainK; trainL; trainM; trainN];
trainLabels = [trainALabels ; trainBLabels; trainCLabels; trainDLabels; trainELabels; trainFLabels; trainGLabels; trainHLabels; trainILabels; trainJLabels; trainKLabels; trainLLabels; trainMLabels; trainNLabels];
testData = [testA ; testB; testC; testD; testE; testF; testG; testH; testI; testJ; testK; testL; testM; testN];
testLabels = [testALabels; testBLabels; testCLabels; testDLabels; testELabels; testFLabels; testGLabels; testHLabels; testILabels; testJLabels; testKLabels; testLLabels; testMLabels; testNLabels];
