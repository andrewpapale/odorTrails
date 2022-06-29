% gettfinalcontrol

t2s = nan(1,length(mod_out));

for i = 1:length(mod_out)
    if mod_out(i).success
        t2s(i) = find(mod_out(i).dnT < 2,1,'first');
    else
        t2s(i) = -20;
    end
end