clear all
close all

%% LOAD

% load 2016-07-21_patches_1_2390um
% load 2016-07-21_patches_2_2442um
% load 2016-07-21_patches_3_2490um
% load 2016-07-21_patches_4_254 0um
% load 2016-07-21_patches_5_2595um
% load 2016-07-21_patches_6_2645um
% load 2016-07-21_patches_7_2695um
% load 2016-07-21_patches_8_2750um
% load 2016-07-21_patches_9_2805um
% load 2016-07-21_patches_10_2860um
% load 2016-07-21_patches_11_2915um

% load 2016-08-26_squares_PH_6

dir1='Presentation';
dir2='2016-07-21';
dir3='24.90';
dir4=strcat(dir2,'\',dir3);
mkdir(dir1,dir4);
dir=strcat(dir1,'\',dir4);

thrX=0.78;
thrY=0.78;

global Mx My MaskX MaskY
Mx=size(x);     My=size(y);
Mx=Mx(2);       My=My(2);

i=1; 

%% Gamma calculation

Vrms=0.255; % Volts

Vpp=2*sqrt(2)*Vrms;
scale=0.2218;
TwoDeltaL=Vpp/scale;
DeltaL=TwoDeltaL/2;
gamma=4*pi*DeltaL/10.6;

J1=besselj(1,gamma);
J2=besselj(2,gamma);

C1=1/J1;
C2=1/J2;

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
suptitle('AFM TOPOGRAPHY')
saveas(F1,strcat(dir,'\','1.Topography.png'))

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

% 1ST sideband
SideB1=sqrt(SideB1x.^2+SideB1y.^2);
SideB1phase=MaskX.*atan2(SideB1y,SideB1x)*180/pi;

% 2ND sideband
SideB2=sqrt(SideB2x.^2+SideB2y.^2);
SideB2phase=MaskX.*atan2(SideB2y,SideB2x)*180/pi;

F2=figure('units','normalized','outerposition',[0 0 1 1]);

% 1ST sideband
% suptitle(strcat(dir3,' (Non Rotated)'))
% suptitle('SIDEBANDS (METHOD 1)')
FigSB1x=subplot(2,3,1);
imagesc(x,y,SideB1x)
title('1st Sideband X')
colormap(FigSB1x,hot)
bar=colorbar
xlabel(bar,'(V)')
axis square
xlabel('x (\mum)')
ylabel('y (\mum)');

FigSB1y=subplot(2,3,2);
imagesc(x,y,SideB1y)
title('1st Sideband Y')
colormap(FigSB1y,hot)
bar=colorbar
xlabel(bar,'(V)')
axis square
xlabel('x (\mum)')
ylabel('y (\mum)');


FigSB1=subplot(2,3,3);
imagesc(x,y,SideB1)
title('1st Sideband Modulus')
colormap(FigSB1,hot)
bar=colorbar
xlabel(bar,'(V)')
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')

% FigSB1ph=subplot(2,4,4);
% imagesc(x,y,SideB1phase)
% title('1st Sideband Phase')
% colormap(FigSB1ph,rwb)
% colorbar
% axis square
% xlabel('x (\mum)')
% ylabel('y (\mum)')

% 2ND sideband
FigSB2x=subplot(2,3,4);
imagesc(x,y,SideB2x)
title('2nd Sideband X')
colormap(FigSB2x,hot)
bar=colorbar
xlabel(bar,'(V)')
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')

FigSB2y=subplot(2,3,5);
imagesc(x,y,SideB2y)
title('2nd Sideband Y')
colormap(FigSB2y,hot)
bar=colorbar
xlabel(bar,'(V)')
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')



FigSB2=subplot(2,3,6);
imagesc(x,y,SideB2)
title('2nd Sideband Modulus')
colormap(FigSB2,hot)
bar=colorbar
xlabel(bar,'(V)')
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')

% FigSB2ph=subplot(2,4,8);
% imagesc(x,y,SideB2phase)
% title('1st Sideband Phase')
% colormap(FigSB2ph,rwb)
% colorbar
% axis square
% xlabel('x (\mum)')
% ylabel('y (\mum)')
saveas(F2,strcat(dir,'\','2.SideBands_Method1.png'))

%% NEAR FIELD (Sideband Modulus)

Tau=C1*SideB2-1i*C2*SideB1;

ModulTau=abs(Tau);
PhaseTau=atan2(imag(Tau),real(Tau)).*MaskX*180/pi;
% PhaseTau=angle(Tau).*MaskX*180/pi;

F3=figure('units','normalized','outerposition',[0 0 1 1]);
% suptitle('NEAR-FIELD (METHOD 1)')
Fig1=subplot(1,2,1);
imagesc(x,y,ModulTau)
title('|E|')
colormap(Fig1,hot)
bar=colorbar
xlabel(bar,'(V)')
axis square
xlabel('x (\mum)')
ylabel('y (\mum)');

Fig2=subplot(1,2,2);
imagesc(x,y,PhaseTau)
title('E Phase')
colormap(Fig2,rwb)
bar=colorbar
xlabel(bar,'(deg)')
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')
saveas(F3,strcat(dir,'\','3.Tau.png'))

%% NEAR FIELD (Only X component)

Tau=C1*SideB2x-1i*C2*SideB1x;

ModulTau=abs(Tau);
PhaseTau=atan2(imag(Tau),real(Tau)).*MaskX*180/pi;
% PhaseTau=angle(Tau).*MaskX*180/pi;

F4=figure('units','normalized','outerposition',[0 0 1 1]);
suptitle('NEAR-FIELD (X COMPONENT)')
Fig1=subplot(1,2,1);
imagesc(x,y,ModulTau)
title('|E|')
colormap(Fig1,hot)
bar=colorbar
xlabel(bar,'(V)')
axis square
xlabel('x (\mum)')
ylabel('y (\mum)');

Fig2=subplot(1,2,2);
imagesc(x,y,PhaseTau)
title('E Phase')
colormap(Fig2,rwb)
bar=colorbar
xlabel(bar,'(deg)')
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')
saveas(F4,strcat(dir,'\','4.Tau.png'))

%% ROTATION

[X1,Y1,X1rot,Y1rot,X1subs,Y1subs,X1subsrot,Y1subsrot,MaxEigenVect1,MinEigenVect1,SideB1xrot,SideB1yrot]=Covariance(K,SideB1x,SideB1y);

[X2,Y2,X2rot,Y2rot,X2subs,Y2subs,X2subsrot,Y2subsrot,MaxEigenVect2,MinEigenVect2,SideB2xrot,SideB2yrot]=Covariance(K,SideB2x,SideB2y);

F5=figure('units','normalized','outerposition',[0 0 1 1]);
subplot(2,2,1)
scatter(X1subs,Y1subs,'.','b')
hold on
scatter(X1,Y1,'.','r')
quiver(0,0,MaxEigenVect1(1),MaxEigenVect1(2),'-m','LineWidth',2);
quiver(0,0,MinEigenVect1(1),MinEigenVect1(2), '-g', 'LineWidth',2);
hold off
axis equal
title('1st Sideband')
legend('substate','patch')
xlabel('X (V)')
ylabel('Y (V)')

subplot(2,2,2)
scatter(X1rot,Y1rot,'.','r')
% hold on
% plot(Xellip1,Yellip1)
% hold off
axis equal
title('1st Sideband (rotated)')
xlabel('X (V)')
ylabel('Y (V)')

subplot(2,2,3)
scatter(X2subs,Y2subs,'.','b')
hold on
scatter(X2,Y2,'.','r')
quiver(0,0,MaxEigenVect2(1),MaxEigenVect2(2),'-m','LineWidth',2);
quiver(0,0,MinEigenVect2(1),MinEigenVect2(2), '-g', 'LineWidth',2);
hold off
axis equal
title('2nd Sideband')
legend('substate','patch')
xlabel('X (V)')
ylabel('X (V)')

subplot(2,2,4)
scatter(X2rot,Y2rot,'.','r')
% hold on
% plot(Xellip2,Yellip2)
% hold off
axis equal
title('2nd Sideband (rotated)')
legend('patch')
xlabel('X (V)')
ylabel('Y (Y)')
% suptitle('SIDEBANDS')
saveas(F5,strcat(dir,'\','5.RotationScatter.png'))

%% SIDEBANDS - ROTATED

% 1ST sideband
SideB1rot=sqrt(SideB1xrot.^2+SideB1yrot.^2);
SideB1rotphase=MaskX.*atan2(SideB1yrot,SideB1xrot)*180/pi;

% 2ND sideband
SideB2rot=sqrt(SideB2xrot.^2+SideB2yrot.^2);
SideB2rotphase=MaskX.*atan2(SideB2yrot,SideB2xrot)*180/pi;

% 1ST sideband

F6=figure('units','normalized','outerposition',[0 0 1 1]);
% suptitle('SIDEBANDS (METHOD 2)')
FigSB1xrot=subplot(2,3,1);
imagesc(x,y,SideB1xrot)
title('1st Sideband X')
colormap(FigSB1xrot,hot)
colorbar
axis square
bar=colorbar
xlabel(bar,'(V)')
xlabel('x (\mum)')
ylabel('y (\mum)');

FigSB1yrot=subplot(2,3,2);
imagesc(x,y,SideB1yrot)
title('1st Sideband Y')
colormap(FigSB1yrot,hot)
bar=colorbar
xlabel(bar,'(V)')
axis square
xlabel('x (\mum)')
ylabel('y (\mum)');

FigSB1rot=subplot(2,3,3);
imagesc(x,y,SideB1rot)
title('1st Sideband Modulus')
colormap(FigSB1rot,hot)
bar=colorbar
xlabel(bar,'(V)')
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')

% FigSB1rotph=subplot(2,4,4);
% imagesc(x,y,SideB1rotphase)
% title('1st Sideband Phase')
% colormap(FigSB1rotph,rwb)
% colorbar
% axis square
% xlabel('x (\mum)')
% ylabel('y (\mum)')

% 2ND sideband

FigSB2xrot=subplot(2,3,4);
imagesc(x,y,SideB2xrot)
title('2nd Sideband X')
colormap(FigSB2xrot,hot)
bar=colorbar
xlabel(bar,'(V)')
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')

FigSB2yrot=subplot(2,3,5);
imagesc(x,y,SideB2yrot)
title('2nd Sideband Y')
colormap(FigSB2yrot,hot)
bar=colorbar
xlabel(bar,'(V)')
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')

FigSB2rot=subplot(2,3,6);
imagesc(x,y,SideB2rot)
title('2nd Sideband Modulus')
colormap(FigSB2rot,hot)
bar=colorbar
xlabel(bar,'(V)')
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')

% FigSB2rotph=subplot(2,4,8);
% imagesc(x,y,SideB2rotphase)
% title('1st Sideband Phase')
% colormap(FigSB2rotph,rwb)
% colorbar
% axis square
% xlabel('x (\mum)')
% ylabel('y (\mum)')
saveas(F6,strcat(dir,'\','6.SideBands_Method2.png'))

%% NEAR FIELD - ROTATED

Taurot=C1*SideB2xrot-1i*C2*SideB1xrot;

ModulTaurot=abs(Taurot);
PhaseTaurot=atan2(imag(Taurot),real(Taurot)).*MaskX*180/pi;
% PhaseTaurot=angle(Taurot).*MaskX*180/pi;

F7=figure('units','normalized','outerposition',[0 0 1 1]);
% suptitle('NEAR-FIELD METHOD 2')
Fig1=subplot(1,2,1);
imagesc(x,y,ModulTaurot)
title('|E| (Rotated)')
colormap(Fig1,hot)
bar=colorbar
xlabel(bar,'(V)')
axis square
xlabel('x (\mum)')
ylabel('y (\mum)');

Fig2=subplot(1,2,2);
imagesc(x,y,PhaseTaurot)
title('E Phase (Rotated)')
colormap(Fig2,rwb)
bar=colorbar
xlabel(bar,'(deg)')
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')
saveas(F7,strcat(dir,'\','7.TauRotated.png'))

%% NEAR FIELD IN COMPLEX PLANE VISUALIZATION - NON ROTATED

[X1,Y1,X1rot,Y1rot,X1subs,Y1subs,X1subsrot,Y1subsrot,MaxEigenVect1,MinEigenVect1,SideB1xrotNF,SideB1yrotNF]=Covariance(K,C2*SideB2x,-C1*SideB1x);

FTau1=figure('units','normalized','outerposition',[0 0 1 1]);
scatter(X1subs,Y1subs,'.','b')
hold on
scatter(X1,Y1,'.','r')
quiver(0,0,MaxEigenVect1(1),MaxEigenVect1(2),'-m','LineWidth',2);
quiver(0,0,MinEigenVect1(1),MinEigenVect1(2), '-g', 'LineWidth',2);
hold off
axis equal
title('NEAR-FIELD (X COMPONENT - METHOD 1)')
legend('substate','patch')
xlabel('S_2_,_2 (V)')
ylabel('S_2_,_1 (V)')
saveas(FTau1,strcat(dir,'\','TauScatter.png'))
%% NEAR FIELD IN COMPLEX PLANE VISUALIZATION - ROTATED

[X1,Y1,X1rot,Y1rot,X1subs,Y1subs,X1subsrot,Y1subsrot,MaxEigenVect1,MinEigenVect1,SideB1xrotNF,SideB1yrotNF]=Covariance(K,C2*SideB2xrot,-C1*SideB1xrot);

FTau2=figure('units','normalized','outerposition',[0 0 1 1]);
scatter(X1subs,Y1subs,'.','b')
hold on
scatter(X1,Y1,'.','r')
quiver(0,0,MaxEigenVect1(1),MaxEigenVect1(2),'-m','LineWidth',2);
quiver(0,0,MinEigenVect1(1),MinEigenVect1(2), '-g', 'LineWidth',2);
hold off
axis equal
title('NEAR-FIELD (X COMPONENT - METHOD 2)')
legend('substate','patch')
xlabel('S_2_,_2 (V)')
ylabel('S_2_,_1 (V)')
saveas(FTau2,strcat(dir,'\','TauScatterRotated.png'))