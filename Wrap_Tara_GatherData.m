%Wrap_Tara_GatherData

Dtop = dir;
nDtop = length(Dtop);

A = {};
iC = 1;
iM = 0;
for iT=3:nDtop-1
    
    cd(Dtop(iT).name);
    
    Dpv = dir('*-PV_cellcounts_19-Nov-2019*');
    Dvva = dir('*-VVA_cellcounts_19-Nov-2019*');
    Dcolo = dir('*colocalization_with_layer*');
    nDpv = length(Dpv);
    nDvva = length(Dvva);
    nDcolo = length(Dcolo);
    
    assert(nDpv==nDvva,'mismatch in pv and vva');
    assert(nDpv==nDcolo,'mismatch in pv and colo');
    
    
    nD = nDpv;
    iM = iM+1;
    
    for iD=1:nD
        
        load(Dpv(iD).name,'cell_pos1');
        nPV = length(cell_pos1);
        load(Dvva(iD).name,'cell_pos1');
        nVVA = length(cell_pos1);
        load(Dcolo(iD).name,'nCoLo','nCoLo14','nCoLo56','nPV14','nPV56','nVVA14','nVVA56');
        
        tempStr = strsplit(Dpv(iD).name,'-');
        
        A{iC,1}=tempStr{1};
        A{iC,2}=tempStr{2};
        A{iC,3}= iM;
        A{iC,4}=nPV;
        A{iC,5}=nPV14;
        A{iC,6}=nPV56;
        A{iC,7}=nVVA;
        A{iC,8}=nVVA14;
        A{iC,9}=nVVA56;
        A{iC,10}=nCoLo;
        A{iC,11}=nCoLo14;
        A{iC,12}=nCoLo56;
        iC = iC+1;
        fprintf('%d/%d \n',iD,nD);
    end
    
    cd('C:\Users\papalea\Documents\Data\Imaging Tara\Re-PV-VVA-withmat');
end
