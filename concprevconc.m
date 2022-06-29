% conc/prevconc
cpc = nan(size(conc0));
for iM=1:nM 
    for iS=1:nS 
        if ~isnan(conc0(iM,iS)) & ~isnan(prevconc(iM,iS)); 
            if conc0(iM,iS)==0 & prevconc(iM,iS)==0
                cpc(iM,iS)=1; 
            elseif conc0(iM,iS)==0 & prevconc(iM,iS)==1 
                cpc(iM,iS)=2;
            elseif conc0(iM,iS)==0 & prevconc(iM,iS)==2
                cpc(iM,iS)=3; 
            elseif conc0(iM,iS)==1 & prevconc(iM,iS)==0
                cpc(iM,iS)=4; 
            elseif conc0(iM,iS)==1 & prevconc(iM,iS)==1
                cpc(iM,iS)=5; 
            elseif conc0(iM,iS)==1 & prevconc(iM,iS)==2
                cpc(iM,iS)=6; 
            elseif conc0(iM,iS)==2 & prevconc(iM,iS)==0
                cpc(iM,iS)=7; 
            elseif conc0(iM,iS)==2 & prevconc(iM,iS)==1
                cpc(iM,iS)=8; 
            elseif conc0(iM,iS)==2 & prevconc(iM,iS)==2
                cpc(iM,iS)=9;
            end
        end
    end
end