% doWelchDline

nM = max(mouse);
nS = max(sess);

winsize = 150;
winstep = 0.85;
winparam = 16;
freqset = linspace(1,10,50);
Fs = 50;

S = [];
S2 = [];
times = [];
times2 = [];
sess1 = [];
mouse1 = [];
sess2 = [];
mouse2 = [];
for iM=1:nM
    for iS=1:nS
        k = mouse==iM & sess==iS;
        
        if sum(k)>0 % GO
            edge0 = edge(k);
            %isstopped0 = isstopped(k);
            var0 = unwrap(polar_phi(k));
            var0(edge0) = nan;
            %var2 = log10dphi(k);
            %var0(edge0 | isstopped0)=nan;
            %var2(edge0 | isstopped0)=nan;
            %t0 = t(k);
            %dt0 = cat(1,nan,diff(t0));
            if ~all(isnan(var0))
                good = find(~isnan(var0)); % should be the same
                %good2 = find(~isnan(var2));
                var0 = interp1(good,var0(good),(1:length(var0))','nearest');
                %var2 = interp1(good2,var2(good2),(1:length(var2))','nearest');
                var0(isnan(var0))=[];
                %var2(isnan(var2))=[];
                winsize0 = min(length(var0),winsize);
                %winsize2 = min(length(var2),winsize);
                winstep0 = round(winstep*winsize0);
                %winstep2 = round(winstep*winsize2);
                
                
                try
                    [~,freq,times0,S0,~,~]=spectrogram(var0,kaiser(winsize0,winparam),winstep0,freqset,Fs);
                    %[~,~,times20,S20,~,~]=spectrogram(var2,kaiser(winsize2,winparam),winstep2,freqset,Fs);
                catch
                    times0 = nan;
                   % times20 = nan;
                    S0 = nan(length(freq),1);
                    %S20 = nan(length(freq),1);
                end
                
                S = cat(2,S,S0);
                %S2 = cat(2,S2,S20);
                times = cat(2,times,times0);
                %times2 = cat(2,times2,times20);
                %if length(times2)~=size(S2,2)
                %    keyboard;
                %end
                sess1 = cat(2,sess1,repmat(iS,[1,length(times0)]));
                mouse1 = cat(2,mouse1,repmat(iM,[1,length(times0)]));
                %sess2 = cat(2,sess2,repmat(iS,[1,length(times20)]));
                %mouse2 = cat(2,mouse2,repmat(iM,[1,length(times20)]));
            end
        end
    end
    disp(iM);
end
            
            
            