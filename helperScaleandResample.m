function resampledFolder = helperScaleandResample(sampleRate,dataFolder)
% This function is only in support of XpwWaveletMLExample. It may change or
% be removed in a future release.
currentDir = pwd;
cd(dataFolder)
a = dir;
a(1:2) = [];
    for j = 1:length(a)
        if(a(j).isdir ==1 )
        cd(dataFolder)
        cd (strcat(pwd,'/', a(j).name))
        sdata = dir('*.mat');
        sinfo = dir('*.info');
            if strcmpi(a(j).name, 'ARRData')
                Freq_sample = 360;
            elseif strcmpi(a(j).name, 'CHFData')
                Freq_sample = 250;
            elseif strcmpi(a(j).name,'NSRData')
                Freq_sample = sampleRate;
            end
            [p,q] = rat(sampleRate/Freq_sample);
            cd(dataFolder)
            cd ../            
            if(j==1)
                mkdir ('ECG_Data');            
            end
            resampledFolder = strcat(pwd,'\','ECG_data');            
            for i=1:length(sdata)
                currentdata = fullfile(dataFolder,a(j).name,sdata(i).name);
                currentinfo = fullfile(dataFolder,a(j).name,sinfo(i).name);
                ecgdata = load(currentdata);                
                fid = fopen(currentinfo,'rt');
                if fid == -1
                    error('Cannot open info file');
                end
                data = textscan(fid,'%d%s%f%f%s','HeaderLines',5,'delimiter','\t');
                base = data{4};
                gain = data{3};
                fclose(fid);
                val = ecgdata.val - base;
                val = val./gain;
                cd(resampledFolder)
                if(i==1)
                mkdir(a(j).name);
                end
                cd (strcat(resampledFolder,'/',a(j).name))    
                resampledData = resample(val(1,:), p,q);                %#ok<NASGU>
                save(strcat(a(j).name,num2str(i),'.mat'),'resampledData')
                resampledData = resample(val(2,:), p,q);                 %#ok<NASGU>
                save(strcat(a(j).name,num2str(length(sdata)+i),'.mat'),'resampledData')                  
            end
        end
    end
cd(currentDir)

end

