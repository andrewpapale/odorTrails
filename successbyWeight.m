% SuccessbyWeight
% 2019-06-10 AndyP
% run getTime2Spot4 first


sfw = [];
imw = [];
wgh = [];
pwgh = [];
nvw = [];
idpw = [];
t2sw = [];
dnidw = [];
iC = 1;
for iM=1:nM
    for iS=1:nS
        k0 = mouse0==iM & sess0==iS;
        if sum(k0(:))>0
            reset0 = trial0(iM,iS)==1;

            if reset0 && iS~=1
                sfw(iC) = nanmean(tempsfw);
                wgh(iC)=weight0(iM,iS);
                pwgh(iC) = pweight0(iM,iS);
                imw(iC)=iM;
                nvw(iC)= nanmean(tempnvw);
                idpw(iC) = nanmean(tempidpw);
                t2sw(iC) = nanmean(tempt2sw(tempt2sw>mint2s));
                dnidw(iC) = nanmean(tempdnidw);
                tempsfw = [];
                tempnvw = [];
                tempidpw = [];
                tempt2sw = [];
                tempdnidw = [];
                iC = iC+1;
            elseif ~reset0 && iS~=1
                tempsfw = cat(1,tempsfw,spotfound(iM,iS));
                tempnvw = cat(1,tempnvw,meannV(iM,iS));
                tempidpw = cat(1,tempidpw,meandphi(iM,iS));
                tempt2sw = cat(1,tempt2sw,t2s(iM,iS));
                tempdnidw = cat(1,tempdnidw,Dnid(iM,iS));
            elseif iS==1
                tempsfw = spotfound(iM,iS);
                tempnvw = meannV(iM,iS);
                tempidpw = meandphi(iM,iS);
                tempt2sw = t2s(iM,iS);
                tempdnidw = Dnid(iM,iS); 
            end
            
        end
    end
end



