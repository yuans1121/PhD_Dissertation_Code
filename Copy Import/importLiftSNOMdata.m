clc
clear all

Folder1='SNOM\2018-09-25\';

N=128;

h=90:10:180;
m=size(h);  m=m(2);
scan='2';

M=zeros(N,6,2*m);

for i=1:m
    hs=num2str(h(i));
    SB1=strcat('lift_',hs,'nm_squares_PH_1SB_',scan);
    SB2=strcat('lift_',hs,'nm_squares_PH_2SB_',scan);
    
    f='.dat';

    File1=strcat(Folder1,SB1,f);
    File2=strcat(Folder1,SB2,f);

    Msb1=zeros(N,6);
    Msb2=zeros(N,6);

    Msb1=dlmread(File1,'\t',10,0);
    Msb2=dlmread(File2,'\t',10,0);

    MSB1(:,:,i)=Msb1;
    MSB2(:,:,i)=Msb2;
end

save 2018-09-25_lift_squares_PH.mat h MSB1 MSB2