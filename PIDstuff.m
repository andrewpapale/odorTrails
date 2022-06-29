NsP = 1:76;
kP = zeros(size(NsP)); 
for iT=1:length(NsP)
    if any(TablePosition==NsP(iT)) 
        kP(iT)=1; 
    end
end


X = [TablePosition,meanV];
X = sort(X,1);
