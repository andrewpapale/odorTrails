xbins = linspace(0,130,130);
%nM = max(mouse);
%nS = max(sess);
mO = nan(length(xbins),3);
dO = nan(length(xbins),3);
k = ~edge & V>0.1 & towards & t2s0;
for iC=0:2
    %for iM=1:nM
    % for iS=1:nS
    k1 = conc==iC;
    if sum(k1)>0
        for iX=1:length(xbins)
            if iX==length(xbins)
                k0 = k & k1 & dnT>xbins(iX);
            else
                k0 = k & k1 & dnT>xbins(iX) & dnT<=xbins(iX+1);
            end
            mO(iX,iC+1)=nanmean(abs(orient(k0)));
            dO(iX,iC+1)=nanstderr(abs(orient(k0)));
        end
        % end
        % end
    end
    disp(iC);
end
