clear all
close all

%% LOAD

load 2016-08-25_squarepatch_4
% load 2016-07-21_patches_2_2442um
% load 2016-07-21_patches_3_2490um
% load 2016-07-21_patches_4_2540um
% load 2016-07-21_patches_5_2595um
% load 2016-07-21_patches_6_2645um
% load 2016-07-21_patches_7_2695um
% load 2016-07-21_patches_8_2750um
% load 2016-07-21_patches_9_2805um
% load 2016-07-21_patches_10_2860um
% load 2016-07-21_patches_11_2915um

dir1='Results';
dir2='2016-08-25';
dir3='4. Patches 5x5';
dir4=strcat(dir2,'\',dir3);
mkdir(dir1,dir4);
dir=strcat(dir1,'\',dir4);

thrX=0.69;
thrY=0.69;

global Mx My MaskX MaskY
Mx=size(x);     My=size(y);
Mx=Mx(2);       My=My(2);

frame=5;
patchsize=16;
Dy=40+2*frame;
Dx=Mx;

i=1; 

%% TOPOGRAPHY 

% X component

TopoX=M(:,:,i);
[MaskX,Kx]=Mask(TopoX,thrX); 

% Y component 

TopoY=M(:,:,i+4);
[MaskY,Ky]=Mask(TopoY,thrY);

if Kx<Ky
    K=Kx;
else
    K=Ky;
end

F1=figure('units','normalized','outerposition',[0 0 1 1]);

% X component

subplot(2,2,1)
imagesc(x,y,TopoX)
title('Topography X')
colormap(gray(256))
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')

subplot(2,2,2)
imagesc(x,y,MaskX)
title('Mask X')
colormap(gray(256))
axis square
xlabel('x (\mum)')
ylabel('y (\mum)');

% Y component 

subplot(2,2,3)
imagesc(x,y,TopoY)
title('Topography Y')
colormap(gray(256))
axis square
xlabel('x (\mum)')
ylabel('y (\mum)');

subplot(2,2,4)
imagesc(x,y,MaskY)
title('Mask Y')
colormap(gray(256))
axis square
xlabel('x (\mum)')
ylabel('y (\mum)');
suptitle(dir3)
saveas(F1,strcat(dir,'\','1.Topography.png'))

% ROW SELECTION

% X component

TopoX=M(:,:,i);
[MaskX,Kx]=Mask(TopoX,thrX); 

RowX=zeros(Dy,Dx,4);

figure('units','normalized','outerposition',[0 0 1 1]);
imagesc(TopoX); title('Select Patches for X Component'); axis square;
[Xxcoord,Xycoord]=ginput;

Xx=round(Xxcoord);
Xy=round(Xycoord);

Xtop=Xy-frame;
Xbottom=Xy+patchsize+frame-1;
Xleft=0;
Xright=Dx;

for j=1:4
    RowX(:,:,j)=M(Xtop:Xbottom,Xleft:Xright,j);
end

[sizeY,sizeX]=size(RowX(:,:,1));
dx=x(2)-x(1);
dy=y(2)-y(1);

xx=0:dx:dx*sizeX;
yy=0:dy:dy*sizeY;

% Y Component

TopoY=M(:,:,i+4);
[MaskY,Ky]=Mask(TopoY,thrY); 

RowY=zeros(Dy,Dx,4);

figure('units','normalized','outerposition',[0 0 1 1]);
imagesc(TopoY); title('Select Patches for Y Component'); axis square;
[Yxcoord,Yycoord]=ginput;

Yx=round(Yxcoord);
Yy=round(Yycoord);

Ytop=Yy-frame;
Ybottom=Yy+patchsize+frame-1;
Yleft=0;
Yright=Dx;

for j=1:4
    RowY(:,:,j)=M(Ytop:Ybottom,Yleft:Yright,j+4);
end

%% ANALYSIS

% 1ST SIDEBAND

SideB1x=RowX(:,:,i+2);  % X component
SideB1y=RowY(:,:,i+2);  % Y component

%2ND SIDEBAND

SideB2x=RowX(:,:,i+3);  % X component
SideB2y=RowY(:,:,i+3);  % Y component

% Phase color map

rm=(0:31)'/32; gm=[rm; 1; flipud(rm)]; rm=[rm; ones(33,1)]; bm=flipud(rm);  % red-white-blue color map for phase
rwb=[rm gm bm];   

%% SIDEBANDS - NON ROTATED

% 1ST sideband
SideB1=sqrt(SideB1x.^2+SideB1y.^2);
SideB1phase=MaskX.*atan2(SideB1y,SideB1x)*180/pi;

% 2ND sideband
SideB2=sqrt(SideB2x.^2+SideB2y.^2);
SideB2phase=MaskX.*atan2(SideB2y,SideB2x)*180/pi;

F2=figure('units','normalized','outerposition',[0 0 1 1]);

% 1ST sideband
suptitle(strcat(dir3,' (Non Rotated)'))
FigSB1x=subplot(2,4,1);
imagesc(x,y,SideB1x)
title('1st Sideband X')
colormap(FigSB1x,hot)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)');

FigSB1y=subplot(2,4,2);
imagesc(x,y,SideB1y)
title('1st Sideband Y')
colormap(FigSB1y,hot)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)');

FigSB1=subplot(2,4,3);
imagesc(x,y,SideB1)
title('1st Sideband Modulus')
colormap(FigSB1,hot)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')

FigSB1ph=subplot(2,4,4);
imagesc(x,y,SideB1phase)
title('1st Sideband Phase')
colormap(FigSB1ph,rwb)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')

% 2ND sideband
FigSB2x=subplot(2,4,5);
imagesc(x,y,SideB2x)
title('2nd Sideband X')
colormap(FigSB2x,hot)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')

FigSB2y=subplot(2,4,6);
imagesc(x,y,SideB2y)
title('2nd Sideband Y')
colormap(FigSB2y,hot)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')

FigSB2=subplot(2,4,7);
imagesc(x,y,SideB2)
title('2nd Sideband Modulus')
colormap(FigSB2,hot)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')

FigSB2ph=subplot(2,4,8);
imagesc(x,y,SideB2phase)
title('1st Sideband Phase')
colormap(FigSB2ph,rwb)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')
saveas(F2,strcat(dir,'\','2.SideBands.png'))

%% NEAR FIELD

Tau=SideB2x-1i*SideB1x;

ModulTau=abs(Tau);
PhaseTau=atan2(-SideB1x,SideB2x).*MaskX*180/pi;

F3=figure('units','normalized','outerposition',[0 0 1 1]);
suptitle(strcat(dir3,' (Non Rotated)'))
Fig1=subplot(1,2,1);
imagesc(x,y,ModulTau)
title('|\tau_n |')
colormap(Fig1,hot)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)');

Fig2=subplot(1,2,2);
imagesc(x,y,PhaseTau)
title('\tau_n Phase')
colormap(Fig2,rwb)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')
saveas(F3,strcat(dir,'\','3.Tau.png'))

%% ROTATION

[X1,Y1,X1rot,Y1rot,X1subs,Y1subs,X1subsrot,Y1subsrot,MaxEigenVect1,MinEigenVect1,SideB1xrot,SideB1yrot]=Covariance(K,SideB1x,SideB1y);

[X2,Y2,X2rot,Y2rot,X2subs,Y2subs,X2subsrot,Y2subsrot,MaxEigenVect2,MinEigenVect2,SideB2xrot,SideB2yrot]=Covariance(K,SideB2x,SideB2y);

F4=figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,2,1)
scatter(X1subs,Y1subs,'.','b')
hold on
scatter(X1,Y1,'.','r')
quiver(0,0,MaxEigenVect1(1),MaxEigenVect1(2),'-m','LineWidth',2);
quiver(0,0,MinEigenVect1(1),MinEigenVect1(2), '-g', 'LineWidth',2);
hold off
axis equal
title('1st Sideband')
xlabel('x')
ylabel('y')

subplot(2,2,2)
scatter(X1rot,Y1rot,'.','r')
% hold on
% plot(Xellip1,Yellip1)
% hold off
axis equal
title('1st Sideband (rotated)')
xlabel('x')
ylabel('y')

subplot(2,2,3)
scatter(X2subs,Y2subs,'.','b')
hold on
scatter(X2,Y2,'.','r')
quiver(0,0,MaxEigenVect2(1),MaxEigenVect2(2),'-m','LineWidth',2);
quiver(0,0,MinEigenVect2(1),MinEigenVect2(2), '-g', 'LineWidth',2);
hold off
axis equal
title('2nd Sideband')
xlabel('x')
ylabel('y')

subplot(2,2,4)
scatter(X2rot,Y2rot,'.','r')
% hold on
% plot(Xellip2,Yellip2)
% hold off
axis equal
title('2nd Sideband (rotated)')
xlabel('x')
ylabel('y')
suptitle(dir3)
saveas(F4,strcat(dir,'\','4.Scatter.png'))

%% SIDEBANDS - ROTATED

% 1ST sideband
SideB1rot=sqrt(SideB1xrot.^2+SideB1yrot.^2);
SideB1rotphase=MaskX.*atan2(SideB1yrot,SideB1xrot)*180/pi;

% 2ND sideband
SideB2rot=sqrt(SideB2xrot.^2+SideB2yrot.^2);
SideB2rotphase=MaskX.*atan2(SideB2yrot,SideB2xrot)*180/pi;

% 1ST sideband

F5=figure('units','normalized','outerposition',[0 0 1 1]);
suptitle(strcat(dir3,' (Rotated)'))
FigSB1xrot=subplot(2,4,1);
imagesc(x,y,SideB1xrot)
title('1st Sideband X')
colormap(FigSB1xrot,hot)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)');

FigSB1yrot=subplot(2,4,2);
imagesc(x,y,SideB1yrot)
title('1st Sideband Y')
colormap(FigSB1yrot,hot)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)');



FigSB1rot=subplot(2,4,3);
imagesc(x,y,SideB1rot)
title('1st Sideband Modulus')
colormap(FigSB1rot,hot)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')

FigSB1rotph=subplot(2,4,4);
imagesc(x,y,SideB1rotphase)
title('1st Sideband Phase')
colormap(FigSB1rotph,rwb)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')

% 2ND sideband

FigSB2xrot=subplot(2,4,5);
imagesc(x,y,SideB2xrot)
title('2nd Sideband X')
colormap(FigSB2xrot,hot)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')

FigSB2yrot=subplot(2,4,6);
imagesc(x,y,SideB2yrot)
title('2nd Sideband Y')
colormap(FigSB2yrot,hot)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')

FigSB2rot=subplot(2,4,7);
imagesc(x,y,SideB2rot)
title('2nd Sideband Modulus')
colormap(FigSB2rot,hot)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')

FigSB2rotph=subplot(2,4,8);
imagesc(x,y,SideB2rotphase)
title('1st Sideband Phase')
colormap(FigSB2rotph,rwb)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')
saveas(F5,strcat(dir,'\','5.SideBandsRotated.png'))

%% NEAR FIELD - ROTATED

Taurot=SideB2xrot-1i*SideB1xrot;

ModulTaurot=abs(Taurot);
PhaseTaurot=atan2(-SideB1xrot,SideB2xrot).*MaskX*180/pi;
% PhaseTaurot=angle(Taurot).*MaskX*180/pi;

F6=figure('units','normalized','outerposition',[0 0 1 1]);
suptitle(strcat(dir3,' (Rotated)'))
Fig1=subplot(1,2,1);
imagesc(x,y,ModulTaurot)
title('|\tau_n | (Rotated)')
colormap(Fig1,hot)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)');

Fig2=subplot(1,2,2);
imagesc(x,y,PhaseTaurot)
title('\tau_n Phase (Rotated)')
colormap(Fig2,rwb)
colorbar
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')
saveas(F6,strcat(dir,'\','6.TauRotated.png'))

%% NEAR FIELD IN COMPLEX PLANE VISUALIZATION - ROTATED

[X1,Y1,X1rot,Y1rot,X1subs,Y1subs,X1subsrot,Y1subsrot,MaxEigenVect1,MinEigenVect1,SideB1xrotNF,SideB1yrotNF]=Covariance(K,SideB2xrot,-SideB1xrot);

FTau=figure();
scatter(X1subs,Y1subs,'.','b')
hold on
scatter(X1,Y1,'.','r')
quiver(0,0,MaxEigenVect1(1),MaxEigenVect1(2),'-m','LineWidth',2);
quiver(0,0,MinEigenVect1(1),MinEigenVect1(2), '-g', 'LineWidth',2);
hold off
axis equal
title('Signal Cloud')
xlabel('S_2_,_2 (V)')
ylabel('S_2_,_1 (V)')
saveas(FTau,strcat(dir,'\','TauScatter_2390.png'))