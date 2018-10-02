clc
clear all

%% SECOND AND THIRD HARMONIC

% Folder1='SNOM\2016-07-26\txt\';
% 
% Position1='squares23H_2_2500_XY_';
% Position2='squares23H_2_2550_XY_';
% Position3='squares23H_2_2600_XY_';
% Position4='squares23H_2_2650_XY_';
% Position5='squares23H_2_2700_XY_';
% Position6='squares23H_2_2750_XY_';
% Position7='squares23H_2_2800_XY_';
% Position8='squares23H_2_2850_XY_';
% Position9='squares23H_2_2900_XY_';
% Position10='squares23H_2_2950_XY_';
% Position11='squares23H_2_3000_XY_';
% Position12='squares23H_2_3050_XY_';
% 
% f='.txt';
% 
% N=128;
% 
% File=strcat(Folder1,Position1,'0',f);
% FileID=fopen(File);
% size=textscan(FileID,'%s','delimiter',' ');
% size=str2double(size{1}{6});
% fclose(FileID);
% 
% x=0:size/(N-1):size;
% 
% y=0:size/(N-1):size;
% 
% M1=zeros(N,N,4);
% M2=zeros(N,N,4);
% M3=zeros(N,N,4);
% M4=zeros(N,N,4);
% M5=zeros(N,N,4);
% M6=zeros(N,N,4);
% M7=zeros(N,N,4);
% M8=zeros(N,N,4);
% M9=zeros(N,N,4);
% M10=zeros(N,N,4);
% M11=zeros(N,N,4);
% M12=zeros(N,N,4);
% 
% for k=0:1:3
%     M1(:,:,k+1)=dlmread(strcat(Folder1,Position1,num2str(k),f),'\t',4,0);
%     M2(:,:,k+1)=dlmread(strcat(Folder1,Position2,num2str(k),f),'\t',4,0);
%     M3(:,:,k+1)=dlmread(strcat(Folder1,Position3,num2str(k),f),'\t',4,0);
%     M4(:,:,k+1)=dlmread(strcat(Folder1,Position4,num2str(k),f),'\t',4,0);
%     M5(:,:,k+1)=dlmread(strcat(Folder1,Position5,num2str(k),f),'\t',4,0);
%     M6(:,:,k+1)=dlmread(strcat(Folder1,Position6,num2str(k),f),'\t',4,0);  
%     M7(:,:,k+1)=dlmread(strcat(Folder1,Position7,num2str(k),f),'\t',4,0);
%     M8(:,:,k+1)=dlmread(strcat(Folder1,Position8,num2str(k),f),'\t',4,0);
%     M9(:,:,k+1)=dlmread(strcat(Folder1,Position9,num2str(k),f),'\t',4,0);
%     M10(:,:,k+1)=dlmread(strcat(Folder1,Position10,num2str(k),f),'\t',4,0);
%     M11(:,:,k+1)=dlmread(strcat(Folder1,Position11,num2str(k),f),'\t',4,0);
%     M12(:,:,k+1)=dlmread(strcat(Folder1,Position12,num2str(k),f),'\t',4,0);
% end
% 
% R1=cat(3,M1(:,:,1),M1(:,:,2),sqrt(M1(:,:,3).^2+M1(:,:,4).^2));
% R2=cat(3,M2(:,:,1),M2(:,:,2),sqrt(M2(:,:,3).^2+M2(:,:,4).^2));
% R3=cat(3,M3(:,:,1),M3(:,:,2),sqrt(M3(:,:,3).^2+M3(:,:,4).^2));
% R4=cat(3,M4(:,:,1),M4(:,:,2),sqrt(M4(:,:,3).^2+M4(:,:,4).^2));
% R5=cat(3,M5(:,:,1),M5(:,:,2),sqrt(M5(:,:,3).^2+M5(:,:,4).^2));
% R6=cat(3,M6(:,:,1),M6(:,:,2),sqrt(M6(:,:,3).^2+M6(:,:,4).^2));
% R7=cat(3,M7(:,:,1),M7(:,:,2),sqrt(M7(:,:,3).^2+M7(:,:,4).^2));
% R8=cat(3,M8(:,:,1),M8(:,:,2),sqrt(M8(:,:,3).^2+M8(:,:,4).^2));
% R9=cat(3,M9(:,:,1),M9(:,:,2),sqrt(M9(:,:,3).^2+M9(:,:,4).^2));
% R10=cat(3,M10(:,:,1),M10(:,:,2),sqrt(M10(:,:,3).^2+M10(:,:,4).^2));
% R11=cat(3,M11(:,:,1),M11(:,:,2),sqrt(M11(:,:,3).^2+M11(:,:,4).^2));
% R12=cat(3,M12(:,:,1),M12(:,:,2),sqrt(M12(:,:,3).^2+M12(:,:,4).^2));
% 
% M=cat(3,M1,M2,M3,M4,M5,M6,M7,M8,M9,M10,M11,M12);
% 
% R=cat(3,R1,R2,R3,R4,R5,R6,R7,R8,R9,R10,R11,R12);
% 
% save 2016-07-26_squares_2and3H.mat M R x y 

%% PSEUDOHETERODYNE MODE

Folder2='SNOM\2018-05-11\txt\';

StructX='discrete_bowtie_LP_150_PH_X_5_';
% StructY='squares_6Y_';

f='.txt';

N=128;

File=strcat(Folder2,StructX,'0',f);
FileID=fopen(File);
size=textscan(FileID,'%s','delimiter',' ');
size=str2double(size{1}{6});
fclose(FileID);

x=0:size/(N-1):size;

y=0:size/(N-1):size;

Mx=zeros(N,N,4);
% My=zeros(N,N,4);

for k=0:1:3
    Mxi=dlmread(strcat(Folder2,StructX,num2str(k),f),'\t',4,0);
    Mx(:,:,k+1)=Mxi;
%     My(:,:,k+1)=dlmread(strcat(Folder2,StructY,num2str(k),f),'\t',4,0);
end

% M=cat(3,Mx,My);
M=Mx;

save 2018-05-11_discrete_bowtie_LP_150_PH_X_5.mat M x y