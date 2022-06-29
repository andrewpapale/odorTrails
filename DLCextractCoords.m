function [frame,blx,bly,blp,brx,bry,brp,flx,fly,flp,frx,fry,frp,bx,by,bp,nx,ny,np,tbx,tby,tbp] = DLCextractCoords(Data);

blx=[];
bly=[];
blp=[];
brx=[];
bry=[];
brp=[];
flx=[];
fly=[];
flp=[];
frx=[];
fry=[];
frp=[];
bx=[];
by=[];
bp=[];
nx=[];
ny=[];
np=[];
tbx=[];
tby=[];
tbp=[];


frame = table2array(Data(:,1));
if size(Data,2)==22
    blx = table2array(Data(:,2));
    bly = table2array(Data(:,3));
    blp = table2array(Data(:,4));
    brx = table2array(Data(:,5));
    bry = table2array(Data(:,6));
    brp = table2array(Data(:,7));
    flx = table2array(Data(:,8));
    fly = table2array(Data(:,9));
    flp = table2array(Data(:,10));
    frx = table2array(Data(:,11));
    fry = table2array(Data(:,12));
    frp = table2array(Data(:,13));
    bx = table2array(Data(:,14));
    by = table2array(Data(:,15));
    bp = table2array(Data(:,16));
    nx = table2array(Data(:,17));
    ny = table2array(Data(:,18));
    np = table2array(Data(:,19));
    tbx = table2array(Data(:,20));
    tby = table2array(Data(:,21));
    tbp = table2array(Data(:,22));
elseif size(Data,2)==10
    bx = table2array(Data(:,2));
    by = table2array(Data(:,3));
    bp = table2array(Data(:,4));
    nx = table2array(Data(:,5));
    ny = table2array(Data(:,6));
    np = table2array(Data(:,7));
    tbx = table2array(Data(:,8));
    tby = table2array(Data(:,9));
    tbp = table2array(Data(:,10));
end