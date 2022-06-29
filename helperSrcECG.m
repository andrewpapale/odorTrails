function ECGData = helperSrcECG(src_dir,signalLength)
% This function is only in support of XpwWaveletMLExample. It may change or
% be removed in a future release.
src_files = create_src(src_dir);
Data = zeros(numel(src_files),signalLength);
for numsig = 1:numel(src_files) 
    datastruct = load(char(src_files{numsig}));
    Data(numsig,:) = datastruct.resampledData(1:signalLength);
end
ECGData.Data = Data;
currentDir = pwd;
cd(src_dir)
cd('ARRdata');
darr = dir([pwd, '\*.mat']);
Narr = numel(darr);
cd(src_dir)
cd('CHFdata');
dchf = dir([pwd, '\*.mat']);
Nchf = numel(dchf);
cd(src_dir)
cd('NSRdata');
dnsr = dir([pwd,'\*.mat']);
Nnsr = numel(dnsr);
cd(src_dir)

Numsubjects = sum([Narr Nchf Nnsr]);
val = (1:Numsubjects)';
Map = containers.Map(val,cellstr([repmat('ARR',Narr,1);...
     repmat('CHF',Nchf,1); repmat('NSR',Nnsr,1)]));

for i =1:length(val)
    Labels(i,:) = Map(val(i)); %#ok<AGROW>
end
ECGData.Labels = cellstr(Labels);

cd(currentDir);
end


function src_files = create_src(directory)
if nargin < 1
    error('Must specify directory!')
end
files = find_files(directory);
if isempty(files)
    error('No data files found in the specified directory!');
end

src_files = files;


end
function files = find_files(directory)
extensions = {'mat'};

dir_list = dir(directory);

files = {};

for k = 1:length(dir_list)
    name = dir_list(k).name;
    
    % Hidden file or current/upper directory? Skip.
    if name(1) == '.'
        continue;
    end
    
    % recurse for all files
    if dir_list(k).isdir
        files = [files find_files(fullfile(directory,name))]; %#ok<AGROW>
    else
        found = 0;
        for l = 1:length(extensions)
            len = length(extensions{l});
            if length(name) > len+1 && ...
                    strcmpi(name(end-len:end),['.' extensions{l}])
                found = l;
                break;
            end
        end
        if found > 0
            files = [files fullfile(directory,name)]; %#ok<AGROW>
        end
    end
end
end
