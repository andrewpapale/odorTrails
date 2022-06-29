% 2019-11-26 AndyP
% Tara_Colocalize in Layer

Dtop = dir;
nDtop = length(Dtop);

for iT=3:nDtop-1
    
    
    
    cd(Dtop(iT).name);
    
    Dc = dir('*colocalization*.mat');
    for iD=1:length(Dc); delete(Dc(iD).name); end
    Dpv = dir('*-PV_cellcounts_19-Nov-2019*');
    Dvva = dir('*-VVA_cellcounts_19-Nov-2019*');
    D14 = dir('*Layer 1-4_CutOut*');
    D56 = dir('*Layer 5-6_CutOut*');
    nDpv = length(Dpv);
    nDvva = length(Dvva);
    nD14 = length(D14);
    nD56 = length(D56);
    
    
    assert(nDpv==nDvva,'mismatch in pv and vva');
    assert(nDpv==nD14,'mismatch in pv and layer 1-4');
    assert(nDpv==nD56,'mismatch in pv and layer 5-6');
    
    nD = nDpv;
    
    
    for iD=1:nD

        load(Dpv(iD).name,'cell_pos1'),
        posPV = cell_pos1;
        load(Dvva(iD).name,'cell_pos1','Radius');
        posVVA = cell_pos1;
        radiusVVA = Radius;
        load(D14(iD).name,'roi');
        roi14 = roi;
        load(D56(iD).name,'roi');
        roi56 = roi;
        
        
        nPV = size(posPV,1);
        nVVA = size(posVVA,1);
        ixCoLo = nan(nPV,nVVA);
        nCoLo = 0;
        nPV14 = 0;
        nVVA14 = 0;
        nCoLo14 = 0;
        nPV56 = 0;
        nVVA56 = 0;
        nCoLo56 = 0;
        both = 0;
        
        [in14PV,on14PV] = inpolygon(posPV(:,1),posPV(:,2),roi14.Position(:,1),roi14.Position(:,2));
        nPV14 = sum(in14PV);
        [in14VVA,on14VVA] = inpolygon(posVVA(:,1),posVVA(:,2),roi14.Position(:,1),roi14.Position(:,2));
        nVVA14 = sum(in14VVA);
        [in56PV,on56PV] = inpolygon(posPV(:,1),posPV(:,2),roi56.Position(:,1),roi56.Position(:,2));
        nPV56 = sum(in56PV);
        [in56VVA,on56VVA] = inpolygon(posVVA(:,1),posVVA(:,2),roi56.Position(:,1),roi56.Position(:,2));
        nVVA56 = sum(in56VVA);
        
        for iVVA=1:nVVA
            ix = find((posPV(:,1)-posVVA(iVVA,1)).^2+(posPV(:,2)-posVVA(iVVA,2)).^2 <= radiusVVA(iVVA).^2);
            if ~isempty(ix)
                ixCoLo(ix,iVVA)=1;
                nCoLo = nCoLo+1;
                if in14VVA(iVVA) && ~in56VVA(iVVA)
                    nCoLo14=nCoLo14+1;
                end
                if in56VVA(iVVA) && ~in14VVA(iVVA)
                    nCoLo56=nCoLo56+1;
                end
                if in14VVA(iVVA) && in56VVA(iVVA)
                    both=both+1;
                end
                %plot(posPV(ix,1),posPV(ix,2),'ms','markersize',10);
                if length(ix)>1
                    warning('multiple cells colocalized with VVA cell %d',iVVA);
                end
            end
        end
        
        dateStr = strsplit(datestr(now),' ');
        tempStr = strsplit(Dpv(iD).name,'-');
        saveStr = strcat(tempStr{1},'-',tempStr{2},'_','colocalization_with_layer','_',dateStr{1},'.mat');
        save(saveStr,'both','nCoLo','ixCoLo','nPV14','nVVA14','nPV56','nVVA56','nCoLo14','nCoLo56');
        
        fprintf('%d/%d \n',iD,nD);
    end
    
    cd('C:\Users\papalea\Documents\Data\Imaging Tara\Re-PV-VVA-withmat');
end