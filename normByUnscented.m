% 2018-10-26 AndyP
% normalize number of trials by nTrials @ C==0 (unscented)
nM = 5;
nS = 105;
nTr = nan(size(sT(:)));

Tr0 = Tr';
sT0 = sT';
C0 = C';
Tr0 = Tr0(:);
sT0 = sT0(:);
C0 = C0(:);
newday = find(Tr0==1);
unscented = find(C0==0);
%U = nan(size(nTr));

for iD=1:length(newday)
    if iD~=length(newday)
        tstop = newday(iD+1)-1;
    else
        tstop = length(nTr);
    end
    nTr(newday(iD):tstop) = sT0(newday(iD):tstop)./(nanmax(sT0));
    %currU = sT0(unscented >= newday(iD) & unscented <= tstop);
%     if length(currU)==1
%         %U(newday(iD):tstop) = repmat(currU,[1,length(newday(iD):tstop)]);
%         
%        % disp(currU);
%     else
%         fprintf('currU>1 or empty %d \n',newday(iD));
%     end
end

%U = U';
nTr = nTr';
nTr = reshape(nTr,[5,105]);
%U = reshape(U,[5, 105]);    
    