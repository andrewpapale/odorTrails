function Wrap_vid2hdf5

DtoDo = dir('*.avi');
Ddone = dir('*.h5');

for iD=1:length(DtoDo)
    Dname = DtoDo(iD).name;
    Dname0 = strsplit(Dname,'.');
    Dname0 = Dname0{1};
    deleteflag = false;
    for iE=1:length(Ddone)
        Ename = Ddone(iE).name;
        Ename = strsplit(Ename,'.');
        Ename = Ename{1};
        if strcmp(Dname0,Ename)
            deleteflag = true;
            delete(Dname);
        break;
        end
    end
    if ~deleteflag % run
        disp(Dname);
        try
            vid2hdf5(Dname);
            delete(Dname);
        catch
            warning('unable to convert video %s',Dname);
        end
    end
end
