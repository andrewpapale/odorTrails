
D = dir('*.avi');
nD = length(D);

for iD=1:nD
    Dname = D(iD).name;
    vid2hdf5_Hooks_chunk(Dname,1);
    vid2hdf5_Hooks_chunk(Dname,7000);
    disp(iD);
end