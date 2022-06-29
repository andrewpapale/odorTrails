spotfound = nan(1,length(mod_out));
for i = 2:length(mod_out)
    if mod_out(i).success
        spotfound(i) = 1;
    else
        spotfound(i) = 0;
    end
end