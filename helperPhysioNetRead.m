function helperPhysioNetRead(dataFolder)
% This function is only intended to support the XpwWaveletMLExample. It
% may change or be removed in a future release.

% Get current directory. Return to current directory after read operation. 
currentDir = pwd;
% Global options
options = weboptions;
options.Timeout = 120;

% Change to input folder. We assume user has write permission and
% sufficient disk space.
cd(dataFolder);
mkdir('rawData');
fprintf('rawData folder created\n');
dataFolder = [dataFolder,'\rawData'];
% Change to rawData folder
cd(dataFolder);
success = mkdir('ARRData');
if success
    cd('ARRData');
    mitdbRead(options);
    fprintf('MIT-BIH Arrhythmia Database Downloaded\n');
elseif ~success
    error('Failed to create ARRData folder');
end
cd(dataFolder);
success = mkdir('CHFData');
if success
    cd('CHFData');
    % 15 records
    chfdbRead(15,options);
    fprintf('BIDMC Congestive Heart Failure Database Downloaded\n');
elseif ~success
    error('Failed to create CHFData folder');
end
cd(dataFolder);
success = mkdir('NSRData');
if success
    cd('NSRData');
    nsrdbRead(options);
    fprintf('MIT-BIH Normal Sinus Rhythm Database Downloaded\n');
elseif ~success
    error('Failed to create NSRData folder');
end
fprintf('All databases downloaded\n');
cd(currentDir);




function chfdbRead(numrecords,options)
for fileInd = 1:numrecords
    fileNumber = sprintf("%02d", fileInd);
    recordName = "chf" + fileNumber;
    fileNameMAT = recordName + "m.mat";
    fileNameINFO = recordName+"m.info";

    webread("http://physionet.org/cgi-bin/atm/ATM?", ...
        "tool", "samples_to_mat", "database", "chfdb", "rbase", recordName, "tdur", "3600", ...
        options);

    url2 = "http://physionet.org/atm/chfdb/" + recordName + "/ecg/0/3600/export/matlab/" + fileNameMAT;
    url3 = "http://physionet.org/atm/chfdb/" + recordName + "/ecg/0/3600/export/matlab/" + fileNameINFO;

     websave(fileNameMAT, url2, options);
     websave(fileNameINFO,url3,options);
end


function mitdbRead(options)
recordName = ["100","101","102","103","104","105","106","107","108","109",...
"111","112","113","114","115","116","117","118","119","121","122","123",...
"124","200","201","202","203","205","207","208","209","210","212","213",...
"214","215","217","219","220","221","222","223","228","230","231","232",...
"233","234"];
for fileInd = 1:numel(recordName)
    fileNameMAT = recordName(fileInd) + "m.mat";
    fileNameINFO = recordName(fileInd)+"m.info";

    webread("http://physionet.org/cgi-bin/atm/ATM?", ...
        "tool", "samples_to_mat", "database", "mitdb", "rbase", recordName(fileInd), "tdur", "3600", ...
        options);

    url2 = "http://physionet.org/atm/mitdb/" + recordName(fileInd) + "/atr/0/3600/export/matlab/" + fileNameMAT;
    url3 = "http://physionet.org/atm/mitdb/" + recordName(fileInd) + "/atr/0/3600/export/matlab/" + fileNameINFO;

     websave(fileNameMAT, url2, options);
     websave(fileNameINFO,url3,options);
end

function nsrdbRead(options)
recordName = ["16265","16272","16273","16420","16483","16539","16773",...
    "16786","16795","17052","17453","18177","18184","19088","19090",...
    "19093","19140","19830"];
for fileInd = 1:numel(recordName)
    fileNameMAT = recordName(fileInd) + "m.mat";
    fileNameINFO = recordName(fileInd)+"m.info";

    webread("http://physionet.org/cgi-bin/atm/ATM?", ...
        "tool", "samples_to_mat", "database", "nsrdb", "rbase", recordName(fileInd), "tdur", "3600", ...
        options);

    url2 = "http://physionet.org/atm/nsrdb/" + recordName(fileInd) + "/atr/0/3600/export/matlab/" + fileNameMAT;
    url3 = "http://physionet.org/atm/nsrdb/" + recordName(fileInd) + "/atr/0/3600/export/matlab/" + fileNameINFO;

     websave(fileNameMAT, url2, options);
     websave(fileNameINFO,url3,options);
end


