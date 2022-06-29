% get lapstart

doTest = false;

nM = max(mouse);
nS = max(sess);
tff = nan(size(f2s));
fff = nan(size(f2s));
fff1 = nan(size(f2s));
tff1 = nan(size(f2s));
xff = nan(size(f2s));
yff = nan(size(f2s));
t2s0 = zeros(size(x));
t2s1 = zeros(size(x));
t2s2 = zeros(size(x));
iC = 1;
for iM=1:nM
    for iS=1:nS
        k = mouse==iM & sess==iS;
        kT = mouseT==iM & sessT==iS;
        
        if sum(k)>0 && sum(kT)>0
            f2sN = f2s(kT);
            f2fN = f2f(kT);
            nTrN = nTr(kT);
            nT = sum(kT);
            tstart = find(k==1,1,'first');
            for iT=1:nT
                tend = f2sN(iT);
                if iT==1
                    tlocstart = 0;
                else
                    tlocstart = f2fN(iT-1);
                end
                if ~isnan(tend)
                    assert(tstart+tend <= find(k==1,1,'last'));
                    k0 = zeros(size(k));
                    k0(tstart+tlocstart:tstart+tend) = 1;
                    x0 = x(k0==1);
                    y0 = y(k0==1);
                    nx0 = nx(k0==1);
                    ny0 = ny(k0==1);
                    t0 = t(k0==1);

                    
                    pp = find((x0-560).^2+(y0-30).^2 < 50.^2,1,'last');
                    pp1 = find((x0-560).^2+(y0-30).^2 > 50.^2,1,'first');
                    if isempty(pp)
                        pp = find((x0-560).^2+(y0-30).^2 > 50.^2,1,'first');
                    end
                    
                    fff(iC) = pp+tlocstart;
                    fff1(iC) = pp1+tlocstart;
                    tff(iC) = t0(pp);
                    tff1(iC) = t0(pp1);
                    xff(iC) = x0(pp);
                    yff(iC) = y0(pp);
                    
%                     assert((f2fN(iT)-fff(iC) > 0) | isnan(f2fN(iT)));
%                     assert((f2sN(iT)-fff(iC) > 0) | isnan(f2sN(iT)));
                    if f2sN(iT)-fff(iC) < 0
                        keyboard;
                    end
                    if f2sN(iT)-fff1(iC) < 0
                        keyboard;
                    end

                    t2s2(tstart+pp1+tlocstart:tstart+tend+1)=repmat(1,[1,length(tstart+pp1+tlocstart:tstart+tend+1)]); %#ok<REPMAT>
                    t2s0(tstart+pp1+tlocstart:tstart+tend+1)= repmat(nTrN(iT),[1,length(tstart+tlocstart+pp1:tstart+tend+1)]);
                    t2s1(tstart+pp+tlocstart:tstart+tend+1)= repmat(1,[1,length(tstart+tlocstart+pp:tstart+tend+1)]); %#ok<REPMAT>
                    if doTest
                        F = figure(1); clf;
                        %plot(x(k),y(k),'k.');
                        hold on;
                        x1 = x(k);
                        y1 = y(k);
                        title(mat2str(iT),'fontsize',24);
                        plot(x(k0==1),y(k0==1),'k.');
                        plot(x0(pp),y0(pp),'rx','markersize',10);
                        plot(x1(f2sN(iT)),y1(f2sN(iT)),'gx','markersize',10);
                        plot(x0(pp:end),y0(pp:end),'b.');
                        set(gca,'YLim',[0 480],'XLim',[0 570]);
                        viscircles([560,30],50);
                        pause;
                    end
                    
                    
                end
                iC = iC+1;
            end
        end
        disp(iS);
    end
    disp(iM);
end