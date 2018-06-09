clear all
close all

%% LOAD

load 2016-07-21_patches_11_2915um

dir1='Results';
dir2='2016-07-21';
dir3='Mean 29.15 test';
dir4=strcat(dir2,'\',dir3);
mkdir(dir1,dir4);
dir=strcat(dir1,'\',dir4);

global Mx My MaskX MaskY
global Dy Dx frame Xtop Xbottom Xleft Xright Ytop Ybottom Yleft Yright
Mx=size(x);     My=size(y);
Mx=Mx(2);       My=My(2);

frame=5;
Dy=40+2*frame;
Dx=Dy;

thrX=0.87;
thrY=0.87;

i=1;

gamma=2.63;
J1=besselj(1,gamma);
J2=besselj(2,gamma);
%% TOPOGRAPHY 

% X component

TopoX=M(:,:,i);
[MaskX,Kx]=Mask(TopoX,thrX); 

PatchTopoX=zeros(Dy,Dx,4);

for j=1:4
    figure('units','normalized','outerposition',[0 0 1 1]);
    imagesc(TopoX); title('Select Patches for X Component'); axis square;
    [Xxcoord(j),Xycoord(j)]=ginput;
end

Xx=round(Xxcoord);
Xy=round(Xycoord);

Xtop=Xy-frame;
Xbottom=Xy+40+frame-1;
Xleft=Xx-frame;
Xright=Xx+40+frame-1;

for j=1:4
    PatchTopoX(:,:,j)=TopoX(Xtop(j):Xbottom(j),Xleft(j):Xright(j));
end

[sizeY,sizeX]=size(PatchTopoX(:,:,1));
dx=x(2)-x(1);
dy=y(2)-y(1);

xx=0:dx:dx*sizeX;
yy=0:dy:dy*sizeY;

for j=1:4
    F1=figure(1);
    imagesc(xx,yy,PatchTopoX(:,:,j))
    axis square;
%     saveas(F1,strcat(dir,'\Patch ',num2str(j),'.png'))
end

% Y Component

TopoY=M(:,:,i+4);
[MaskY,Ky]=Mask(TopoY,thrY); 

PatchTopoY=zeros(Dy,Dx,4);

for j=1:4
    figure('units','normalized','outerposition',[0 0 1 1]);
    imagesc(TopoY); title('Select Patches for Y Component'); axis square;
    [Yxcoord(j),Yycoord(j)]=ginput;
end

Yx=round(Yxcoord);
Yy=round(Yycoord);

Ytop=Yy-frame;
Ybottom=Yy+40+frame-1;
Yleft=Yx-frame;
Yright=Yx+40+frame-1;

for j=1:4
    PatchTopoY(:,:,j)=TopoY(Ytop(j):Ybottom(j),Yleft(j):Yright(j));
    F1=figure(1);
    imagesc(xx,yy,PatchTopoY(:,:,j))
    axis square;
%     saveas(F1,strcat(dir,'\Patch ',num2str(j),'.png'))
end

close all

F1=figure('units','normalized','outerposition',[0 0 1 1]);
subplot(1,2,1)
imagesc(x,y,TopoX)
title('Topography X')
colormap(gray(256))
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')

subplot(1,2,2)
imagesc(x,y,TopoY)
title('Topography Y')
colormap(gray(256))
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')
% saveas(F1,strcat(dir,'\','1.Topography.png'))

%% ANALYSIS

% 1ST SIDEBAND

SideB1x=M(:,:,i+2);  % X component
SideB1y=M(:,:,i+6);  % Y component

%2ND SIDEBAND

SideB2x=M(:,:,i+3);  % X component
SideB2y=M(:,:,i+7);  % Y component

% Phase color map

rm=(0:31)'/32; gm=[rm; 1; flipud(rm)]; rm=[rm; ones(33,1)]; bm=flipud(rm);  % red-white-blue color map for phase
rwb=[rm gm bm];

%% SIDEBANDS - NON ROTATED

% 1ST Sideband

F2=figure('units','normalized','outerposition',[0 0 1 1]);
% suptitle(strcat(dir3,' : 2nd Harmonic'))
FigSB1x=subplot(2,2,1);
imagesc(x,y,SideB1x)
title('1st Sideband X')
colormap(FigSB1x,hot)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)');

FigSB1y=subplot(2,2,2);
imagesc(x,y,SideB1y)
title('1st Sideband Y')
colormap(FigSB1y,hot)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)');

% 2ND sideband

FigSB2x=subplot(2,2,3);
imagesc(x,y,SideB2x)
title('2nd Sideband X')
colormap(FigSB2x,hot)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')

FigSB2y=subplot(2,2,4);
imagesc(x,y,SideB2y)
title('2nd Sideband Y')
colormap(FigSB2y,hot)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')
saveas(F2,strcat(dir,'\','2.SideBands.png'))

%% MEAN SIDEBANDS - NON ROTATED

% 1ST Sideband X

[SideB1xMean,SideB1yMean,SideB1xStd,SideB1yStd,SideB1xCV,SideB1yCV,MeanSideB1xyStd]=MeanAndStd(SideB1x,SideB1y);

F3=figure('units','normalized','outerposition',[0 0 1 1]);
% suptitle(strcat(dir3,'\mum : 2nd Harmonic (1st Sideband) Mean'))  
FigSB1xMean=subplot(2,3,1);
imagesc(xx,yy,SideB1xMean)
title('1st Sideband X Mean')
colormap(FigSB1xMean,hot)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)');

FigSB1xStd=subplot(2,3,2);
imagesc(xx,yy,SideB1xStd)
title('1st Sideband X Std')
colormap(FigSB1xStd,bone)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)');

FigSB1xCV=subplot(2,3,3);
imagesc(xx,yy,log10(SideB1xCV))
title('1st Sideband X log(Coeff of Variation)')
colormap(FigSB1xCV,bone)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)');

% 1ST Sideband Y

FigSB1yMean=subplot(2,3,4);
imagesc(xx,yy,SideB1yMean)
title('1st Sideband Y Mean')
colormap(FigSB1yMean,hot)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')

FigSB1yStd=subplot(2,3,5);
imagesc(xx,yy,SideB1yStd)
title('1st Sideband Y Std')
colormap(FigSB1yStd,bone)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')

FigSB1yCV=subplot(2,3,6);
imagesc(xx,yy,log10(SideB1yCV))
title('1st Sideband Y log(Coeff of Variation)')
colormap(FigSB1yCV,bone)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)');
saveas(F3,strcat(dir,'\','3.1stSideBandMean.png'))

% 2ND Sideband X

[SideB2xMean,SideB2yMean,SideB2xStd,SideB2yStd,SideB2xCV,SideB2yCV,MeanSideB2xyStd]=MeanAndStd(SideB2x,SideB2y);

F4=figure('units','normalized','outerposition',[0 0 1 1]);
% suptitle(strcat(dir3,'\mum : 2nd Harmonic (2nd Sideband) Mean'))
FigSB2xMean=subplot(2,3,1);
imagesc(xx,yy,SideB2xMean)
title('2nd Sideband X Mean')
colormap(FigSB2xMean,hot)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)');

FigSB2xStd=subplot(2,3,2);
imagesc(xx,yy,SideB2xStd)
title('2nd Sideband X Std')
colormap(FigSB2xStd,bone)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)');

FigSB2xCV=subplot(2,3,3);
imagesc(xx,yy,log10(SideB2xCV))
title('2nd Sideband X log(Coeff of Variation)')
colormap(FigSB2xCV,bone)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)');

% 2DN Sideband Y

FigSB2yMean=subplot(2,3,4);
imagesc(xx,yy,SideB2yMean)
title('2nd Sideband Y Mean')
colormap(FigSB2yMean,hot)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')

FigSB2yStd=subplot(2,3,5);
imagesc(xx,yy,SideB2yStd)
title('2nd Sideband Y Std')
colormap(FigSB2yStd,bone)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')

FigSB2yCV=subplot(2,3,6);
imagesc(xx,yy,log10(SideB2yCV))
title('2nd Sideband Y log(Coeff of Variation)')
colormap(FigSB2yCV,bone)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)');
saveas(F4,strcat(dir,'\','3.2ndSideBandMean.png'))

MeanSideB1xyStd
MeanSideB2xyStd

%------------------------------

SideB1=sqrt(SideB1x.^2+SideB1y.^2);
SideB1phase=MaskX.*atan2(SideB1y,SideB1x)*180/pi;

SideB2=sqrt(SideB2x.^2+SideB2y.^2);
SideB2phase=MaskX.*atan2(SideB2y,SideB2x)*180/pi;

%% NEAR FIELD - NON ROTATED (modulus then mean)
% First calculate the modulus then calculate the mean

Tau=SideB2x/J2+1i*SideB1x/J1;

ModulTau=abs(Tau);
% PhaseTau=atan2(SideB1,SideB2).*MaskX*180/pi;
PhaseTau=atan2(SideB1x/J1,SideB2x/J2).*MaskX*180/pi;

Xtop=round(mean([Xtop;Ytop]));
Xbottom=round(mean([Xbottom;Ybottom]));
Xleft=round(mean([Xleft;Yleft]));
Xright=round(mean([Xright;Yright]));

Ytop=Xtop;
Ybottom=Xbottom;
Yleft=Xleft;
Yright=Xright;

[ModulTauMean,PhaseTauMean,ModulTauStd,PhaseTauStd,ModulTauCV,PhaseTauCV,ModulPhaseMeanStd]=MeanAndStd(ModulTau,PhaseTau);

F5=figure('units','normalized','outerposition',[0 0 1 1]);
% suptitle(strcat(dir3,'\mum :\tau Modulus and Phase Mean'))

% Modulus
Fig1=subplot(2,3,1);
imagesc(xx,yy,ModulTauMean)
title('Mean |\tau_n |')
colormap(Fig1,hot)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')

Fig2=subplot(2,3,2);
imagesc(xx,yy,ModulTauStd)
title('Std |\tau_n |')
colormap(Fig2,bone)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')

Fig3=subplot(2,3,3);
imagesc(xx,yy,log10(ModulTauCV))
title('|\tau_n | log(Coeff of Variation)')
colormap(Fig3,bone)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')

% Phase
Fig4=subplot(2,3,4);
imagesc(xx,yy,PhaseTauMean)
title('Mean \tau_n Phase')
colormap(Fig4,rwb)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')

Fig5=subplot(2,3,5);
imagesc(xx,yy,PhaseTauStd)
title('Std \tau_n Phase')
colormap(Fig5,bone)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')

Fig6=subplot(2,3,6);
imagesc(xx,yy,PhaseTauCV)
title('\tau_n Phase Coeff of Variation')
colormap(Fig6,bone)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')
saveas(F5,strcat(dir,'\','4.ModulusAndPhaseMean.png'))

ModulPhaseMeanStd

%% ROTATION

if Kx<Ky
    K=Kx;
else
    K=Ky;
end

[X1,Y1,X1rot,Y1rot,X1subs,Y1subs,X1subsrot,Y1subsrot,MaxEigenVect1,MinEigenVect1,SideB1xrot,SideB1yrot]=Covariance(K,SideB1x,SideB1y);

[X2,Y2,X2rot,Y2rot,X2subs,Y2subs,X2subsrot,Y2subsrot,MaxEigenVect2,MinEigenVect2,SideB2xrot,SideB2yrot]=Covariance(K,SideB2x,SideB2y);

%% SIDEBANDS - ROTATED

% 1ST sideband

F6=figure('units','normalized','outerposition',[0 0 1 1]);
% suptitle(strcat(dir3,'\mum : 2nd Harmonic (Rotated)'))
FigSB1xrot=subplot(2,2,1);
imagesc(x,y,SideB1xrot)
title('1st Sideband X')
colormap(FigSB1xrot,hot)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)');

FigSB1yrot=subplot(2,2,2);
imagesc(x,y,SideB1yrot)
title('1st Sideband Y')
colormap(FigSB1yrot,hot)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)');

% 2ND sideband

FigSB2xrot=subplot(2,2,3);
imagesc(x,y,SideB2xrot)
title('2nd Sideband X')
colormap(FigSB2xrot,hot)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')

FigSB2yrot=subplot(2,2,4);
imagesc(x,y,SideB2yrot)
title('2nd Sideband Y')
colormap(FigSB2yrot,hot)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')
saveas(F6,strcat(dir,'\','5.SideBandsRotated.png'))

%% MEAN SIDEBANDS - ROTATED

% 1ST Sideband X

[SideB1xrotMean,SideB1yrotMean,SideB1xrotStd,SideB1yrotStd,SideB1xrotCV,SideB1yrotCV,MeanSideB1xyrotStd]=MeanAndStd(SideB1xrot,SideB1yrot);

F7=figure('units','normalized','outerposition',[0 0 1 1]);
% suptitle(strcat(dir3,'\mum : 2nd Harmonic (1st Sideband Rotated) Mean'))  
FigSB1xMean=subplot(2,3,1);
imagesc(xx,yy,SideB1xrotMean)
title('1st Sideband X Mean')
colormap(FigSB1xMean,hot)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)');

FigSB1xStd=subplot(2,3,2);
imagesc(xx,yy,SideB1xrotStd)
title('1st Sideband X Std')
colormap(FigSB1xStd,bone)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)');

FigSB1xCV=subplot(2,3,3);
imagesc(xx,yy,log10(SideB1xrotCV))
title('1st Sideband X log(Coeff of Variation)')
colormap(FigSB1xCV,bone)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)');

% 1ST Sideband Y

FigSB1yMean=subplot(2,3,4);
imagesc(xx,yy,SideB1yrotMean)
title('1st Sideband Y Mean')
colormap(FigSB1yMean,hot)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')

FigSB1yStd=subplot(2,3,5);
imagesc(xx,yy,SideB1yrotStd)
title('1st Sideband Y Std')
colormap(FigSB1yStd,bone)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')

FigSB1yCV=subplot(2,3,6);
imagesc(xx,yy,log10(SideB1yrotCV))
title('1st Sideband Y log(Coeff of Variation)')
colormap(FigSB1yCV,bone)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)');
saveas(F7,strcat(dir,'\','6.1stSideBandRotatedMean.png'))

% 2ND Sideband X

[SideB2xrotMean,SideB2yrotMean,SideB2xrotStd,SideB2yrotStd,SideB2xrotCV,SideB2yrotCV,MeanSideB2xyrotStd]=MeanAndStd(SideB2xrot,SideB2yrot);

F8=figure('units','normalized','outerposition',[0 0 1 1]);
% suptitle(strcat(dir3,'\mum : 2nd Harmonic (2nd Sideband Rotated) Mean'))
FigSB2xMean=subplot(2,3,1);
imagesc(xx,yy,SideB2xrotMean)
title('2nd Sideband X Mean')
colormap(FigSB2xMean,hot)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)');

FigSB2xStd=subplot(2,3,2);
imagesc(xx,yy,SideB2xrotStd)
title('2nd Sideband X Std')
colormap(FigSB2xStd,bone)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)');

FigSB2xCV=subplot(2,3,3);
imagesc(xx,yy,log10(SideB2xrotCV))
title('2nd Sideband X log(Coeff of Variation)')
colormap(FigSB2xCV,bone)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)');

% 2DN Sideband Y

FigSB2yMean=subplot(2,3,4);
imagesc(xx,yy,SideB2yrotMean)
title('2nd Sideband Y Mean')
colormap(FigSB2yMean,hot)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')

FigSB2yStd=subplot(2,3,5);
imagesc(xx,yy,SideB2yrotStd)
title('2nd Sideband Y Std')
colormap(FigSB2yStd,bone)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')

FigSB2yCV=subplot(2,3,6);
imagesc(xx,yy,log10(SideB2yrotCV))
title('2nd Sideband Y log(Coeff of Variation)')
colormap(FigSB2yCV,bone)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)');
saveas(F8,strcat(dir,'\','6.2ndSideBandRotatedMean.png'))

MeanSideB1xyrotStd
SideB1xyPercent=round(100*(MeanSideB1xyrotStd-MeanSideB1xyStd)./MeanSideB1xyStd)

MeanSideB2xyrotStd
SideB2xyPercent=round(100*(MeanSideB2xyrotStd-MeanSideB2xyStd)./MeanSideB2xyStd)

%------------------------------

SideB1rot=sqrt(SideB1xrot.^2+SideB1yrot.^2);
SideB1rotphase=MaskX.*atan2(SideB1yrot,SideB1xrot)*180/pi;

SideB2rot=sqrt(SideB2xrot.^2+SideB2yrot.^2);
SideB2rotphase=MaskX.*atan2(SideB2yrot,SideB2xrot)*180/pi;

%% NEAR FIELD - ROTATED

Taurot=SideB2xrot/J2+1i*SideB1xrot/J1;

ModulTaurot=abs(Taurot);
PhaseTaurot=atan2(SideB1xrot/J1,SideB2xrot/J2).*MaskX*180/pi;

[ModulTaurotMean,PhaseTaurotMean,ModulTaurotStd,PhaseTaurotStd,ModulTaurotCV,PhaseTaurotCV,ModulPhaserotMeanStd]=MeanAndStd(ModulTaurot,PhaseTaurot);

F9=figure('units','normalized','outerposition',[0 0 1 1]);
% suptitle(strcat(dir3,'\mum :\tau Modulus and Phase (Rotated) Mean'))

% Modulus
Fig1=subplot(2,3,1);
imagesc(xx,yy,ModulTaurotMean)
title('Mean |\tau_n |')
colormap(Fig1,hot)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')

Fig2=subplot(2,3,2);
imagesc(xx,yy,ModulTaurotStd)
title('Std |\tau_n |')
colormap(Fig2,bone)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')

Fig3=subplot(2,3,3);
imagesc(xx,yy,log10(ModulTaurotCV))
title('|\tau_n | log(Coeff of Variation)')
colormap(Fig3,bone)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')

% Phase
Fig4=subplot(2,3,4);
imagesc(xx,yy,PhaseTaurotMean)
title('Mean \tau_n Phase')
colormap(Fig4,rwb)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')

Fig5=subplot(2,3,5);
imagesc(xx,yy,PhaseTaurotStd)
title('Std \tau_n Phase')
colormap(Fig5,bone)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')

Fig6=subplot(2,3,6);
imagesc(xx,yy,PhaseTaurotCV)
title('\tau_n Phase Coeff of Variation')
colormap(Fig6,bone)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')
saveas(F9,strcat(dir,'\','7.ModulusAndPhaseRotatedMean.png'))

ModulPhaserotMeanStd

ModulPhasePercent=round(100*(ModulPhaserotMeanStd-ModulPhaseMeanStd)./ModulPhaseMeanStd)

%% Calculating first the mean and then the modulus

Xtop=Xtop(1);
Xbottom=Xbottom(1);
Xleft=Xleft(1);
Xright=Xright(1);

Ytop=Xtop;
Ybottom=Xbottom;
Yleft=Xleft;
Yright=Xright;

MaskX=MaskX(Xtop:Xbottom,Xleft:Xright);

%% NEAR FIELD - NON ROTATED (mean then modulus)
% First calculate the mean then calculate the modulus

Tau2=SideB2xMean/J2+1i*SideB1xMean/J1;

ModulTau2=abs(Tau2);
PhaseTau2=atan2(SideB1xMean/J1,SideB2xMean/J2).*MaskX*180/pi;

% Comparison

ModulTauComp=100*abs(ModulTauMean-ModulTau2)./ModulTauMean;
PhaseTauComp=100*abs((PhaseTauMean-PhaseTau2)./PhaseTauMean);

F10=figure('units','normalized','outerposition',[0 0 1 1]);
% suptitle(strcat(dir3,'\mum :\tau Modulus and Phase Mean'))

% Modulus
FIG1=subplot(2,2,1);
imagesc(xx,yy,ModulTau2)
title('Mean |\tau_n |')
colormap(FIG1,hot)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')

FIG2=subplot(2,2,2);
imagesc(xx,yy,ModulTauComp)
title('Comparison |\tau_n |')
colormap(FIG2,bone)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')

% Phase
FIG3=subplot(2,2,3);
imagesc(xx,yy,PhaseTau2)
title('Mean \tau_n Phase')
colormap(FIG3,rwb)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')

FIG4=subplot(2,2,4);
imagesc(xx,yy,log10(PhaseTauComp))
title('Comparison \tau_n Phase (log)')
colormap(FIG4,bone)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')
saveas(F10,strcat(dir,'\','10.ModulusPhaseMethod2.png'))

%% NEAR FIELD - ROTATED (mean then modulus)
% First calculate the mean then calculate the modulus

Tau2rot=SideB2xrotMean/J2+1i*SideB1xrotMean/J1;

ModulTau2rot=abs(Tau2rot);
PhaseTau2rot=atan2(SideB1xrotMean/J1,SideB2xrotMean/J2).*MaskX*180/pi;

% Comparison

ModulTaurotComp=100*abs(ModulTaurotMean-ModulTau2rot)./ModulTaurotMean;
PhaseTaurotComp=100*abs((PhaseTaurotMean-PhaseTau2rot)./PhaseTaurotMean);

F11=figure('units','normalized','outerposition',[0 0 1 1]);
% suptitle(strcat(dir3,'\mum :\tau Modulus and Phase Mean'))

% Modulus
Fig1rot=subplot(2,2,1);
imagesc(xx,yy,ModulTau2rot)
title('Mean |\tau_n |')
colormap(Fig1rot,hot)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')

Fig2rot=subplot(2,2,2);
imagesc(xx,yy,ModulTaurotComp)
title('Comparison |\tau_n |')
colormap(Fig2rot,bone)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')

% Phase
Fig3rot=subplot(2,2,3);
imagesc(xx,yy,PhaseTau2rot)
title('Mean \tau_n Phase')
colormap(Fig3rot,rwb)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')

Fig4rot=subplot(2,2,4);
imagesc(xx,yy,log10(PhaseTaurotComp))
title('Comparison \tau_n Phase (log)')
colormap(Fig4rot,bone)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')
saveas(F11,strcat(dir,'\','11.ModulusPhaseMethod2Rotated.png'))
