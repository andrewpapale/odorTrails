% Wrap_Tara_CellCounting

Dtop = dir;
nDtop = length(Dtop);


for iT=3:nDtop-1
    
    cd(Dtop(iT).name);
    
    D = dir('*.tif');
    nD = length(D);
    
    for iD=1:nD
        str = strsplit(D(iD).name,'-');
        if strcmp(str{3},'PV_0.tif')  % repeat for PV_1.tif
            counttype = 1;
            skip = 0;
        elseif strcmp(str{3},'VVA_0.tif')
            counttype = 2;
            skip = 0;
        else
            skip = 1;
            counttype=0;
        end
        
        if ~skip
            temp = imread(D(iD).name);
            tempF = localMedSub(temp);
            
            if counttype==1
                temp = double(tempF);
                temp(temp==0)=nan; % 2019-11-05 AndyP
                thrVal0 = nanmean(temp(:))+1.5*nanstd(temp(:));
                Ithr = temp > thrVal0;
                Ithr = medfilt2(Ithr,[3,3]);
                props = {'Area', 'PixelIdxList', 'Centroid', 'Eccentricity', 'Perimeter', 'Solidity'};
                BW = regionprops(Ithr,props);
                Cells = nan(length(BW),1);
                for iP=1:length(BW)
                    if BW(iP).Area < 30
                        Cells(iP)=0;
                    else
                        Cells(iP)=1;
                    end
                end
                Area = [];
                Centroid = [];
                Eccentricity = [];
                Solidity = [];
                for iP=1:length(BW)
                    if Cells(iP)
                        Area = cat(1,Area,BW(iP).Area);
                        Centroid = cat(1,Centroid,BW(iP).Centroid);
                        Eccentricity = cat(1,Eccentricity, BW(iP).Eccentricity);
                        Solidity = cat(1,Solidity, BW(iP).Solidity);
                    end
                end
                
                areas1 = Area;
                cell_pos1 = Centroid;
                eccentricity1 = Eccentricity;
                solidity1 = Solidity;
                
                dateStr = strsplit(datestr(now),' ');
                tempStr = strsplit(D(iD).name,'_');
                saveStr = strcat(tempStr{1},'_','cellcounts','_',dateStr{1},'.mat');
                save(saveStr,'areas1','cell_pos1','eccentricity1','solidity1');
                
                
                
            elseif counttype==2
                temp = double(tempF);
                thrVal0 = nanmean(temp(:))+2.25*nanstd(temp(:));
                temp = temp > thrVal0;
                temp = medfilt2(temp,[3,3]);
                [cell_pos1,Radius] = imfindcircles(temp,[7 10],'Sensitivity',0.9,'EdgeThreshold',0.1,'Method','TwoStage');
                
                dateStr = strsplit(datestr(now),' ');
                tempStr = strsplit(D(iD).name,'_');
                saveStr = strcat(tempStr{1},'_','cellcounts','_',dateStr{1},'.mat');
                save(saveStr,'cell_pos1','Radius');
                
            end
        end
        fprintf('%d/%d \n',iD,nD);
    end
    
    
    cd('C:\Users\papalea\Documents\Data\Imaging Tara\Re-PV-VVA-withmat');
end
