% get orientation from spot

orient = nan(size(x));
borient = nan(size(x));
dO = nan(size(x));
% IF = nan(size(x));
% IA = nan(size(x));
% IP = nan(size(x));

mouse2 = [];
S = [];
freqs = linspace(0.1,10,100);
params.tapers = [3 9];
params.Fs = 50;
%params.fpass = linspace(0.01,10,250);
mouse3 = [];
S1 = [];
nM = max(mouse);
nS = max(sess);
miss = 0;
spotfound = [];
Oezt = [];
Azdphi = [];

for iM=1:nM
    for iS=1:nS
        k = mouse==iM & sess==iS;
        
        if sum(k)>0
            
            dnT0 = dnT(k);
            ezt = find(dnT0<=5,1,'first');
            
            if ~isempty(ezt)
                
                t1 = (ezt-400);
                t2 = (ezt+400);
                
                
                
                zdphi0 = zdphi1(k);
                
                zdphi0 = zdphi0(max(1,t1):min(sum(k),t2));
                
                if t1(1)<0
                    pad = nan(800-length(zdphi0),1);
                    S1 = cat(1,pad,zdphi0);
                elseif t1(end)>length(k)
                    pad = nan(800-length(zdphi0),1);
                    S1 = cat(1,zdphi0,pad);
                end
                
                
                
                Azdphi = cat(2,Azdphi,S1);
                mouse3 = cat(1,mouse3,iM);
                if t2s(iM,iS)>0
                    spotfound = cat(1,spotfound,1);
                elseif t2s(iM,iS)==-20
                    spotfound = cat(1,spotfound,0);
                elseif t2s(iM,iS)==-40
                    spotfound = cat(1,spotfound,-1);
                end
            else
                miss = miss+1;
                disp('did not cross boundary');
            end
            %[IF(k), IA(k), IP(k)]=InstFreq(orient(k));
        end
    end
    disp(iM);
end