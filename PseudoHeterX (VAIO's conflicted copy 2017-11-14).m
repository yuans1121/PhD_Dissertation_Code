clear all
close all

lambda=10.6; %microns

%% LOAD

load 2017-10-30_bowtie_1_p-pol_PH_X_4

dir1='Presentation';
dir2='2017-10-30';
dir3='Bowtie-Scan4';
dir4=strcat(dir2,'\',dir3);
mkdir(dir1,dir4);
dir=strcat(dir1,'\',dir4);

global Mx My MaskX MaskY
Mx=size(x);     My=size(y);
Mx=Mx(2);       My=My(2);

i=1; 



%% TOPOGRAPHY 

% X component

TopoX=M(:,:,i);
N=size(TopoX);
N=N(1);

MinTopo=min(min(TopoX));
MaxTopo=max(max(TopoX));
thrX=0.10*(MaxTopo-MinTopo)+MinTopo;

StrucInd=find(TopoX>thrX);
SubstInd=find(TopoX<thrX);

% [MaskX,Kx]=Mask(TopoX,thrX); 

MaskX=zeros(N,N);
MaskX(StrucInd)=1;
MaskX(SubstInd)=0;

Square=ones(3);
MaskX=conv2(MaskX,Square,'same');

StrucInd=find(MaskX>0);
SubstInd=find(MaskX==0);
MaskX(StrucInd)=1;

F1=figure('units','normalized','outerposition',[0 0 1 1]);

% X component

subplot(1,2,1);
imagesc(x,y,TopoX)
title('Topography X')
colormap(gray(256))
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')

subplot(1,2,2)
imagesc(x,y,MaskX)
title('Mask X')
colormap(gray(256))
axis square
xlabel('x (\mum)')
ylabel('y (\mum)');
saveas(F1,strcat(dir,'\','1.Topography.png'))

%% ANALYSIS

% 1ST SIDEBAND
SideB1x=M(:,:,i+2);  % X component

%2ND SIDEBAND
SideB2x=M(:,:,i+3);  % X component

% Phase color map

rm=(0:31)'/32; gm=[rm; 1; flipud(rm)]; rm=[rm; ones(33,1)]; bm=flipud(rm);  % red-white-blue color map for phase
PhaseColormap=[rm gm bm];   

%% SIDEBANDS

% 1ST sideband
SideB1=SideB1x;

% 2ND sideband
SideB2=SideB2x;

F2=figure('units','normalized','outerposition',[0 0 1 1]);
% suptitle(strcat(dir3,' (Non Rotated)'))
% suptitle('SIDEBANDS (METHOD 1)')

FigSB1=subplot(1,2,1);
imagesc(x,y,SideB1)
title('1st Sideband')
colormap(FigSB1,hot)
bar=colorbar;
xlabel(bar,'(V)')
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')

% 2ND sideband

FigSB2=subplot(1,2,2);
imagesc(x,y,SideB2)
title('2nd Sideband')
colormap(FigSB2,hot)
bar=colorbar;
xlabel(bar,'(V)')
axis square
xlabel('x (\mum)')
ylabel('y (\mum)')
saveas(F2,strcat(dir,'\','2.SideBands.png'))

%% NEAR FIELD

% Gamma calculation

Vrms=0.178; % Volts

V=sqrt(2)*Vrms;
scale=0.2218;
DeltaL=V/scale;
gamma=4*pi*DeltaL/lambda;

J1=besselj(1,gamma);
J2=besselj(2,gamma);

C1=1/J1;
C2=1/J2;

Vrms=(2.63*lambda/(4*pi))*scale/(sqrt(2))

Tau=C2*SideB2-1i*C1*SideB1;

ModulTau=abs(Tau);
PhaseTau=atan2(-C1*SideB1,C2*SideB2);
% PhaseTau=atan2(imag(Tau),real(Tau));
PhaseTauMask=PhaseTau.*MaskX;

F3=figure('units','normalized','outerposition',[0 0 1 1]);
% suptitle('NEAR-FIELD (METHOD 1)')
Fig1=subplot(1,2,1);
imagesc(x,y,ModulTau)
title('|E|')
colormap(Fig1,parula)
bar=colorbar;
xlabel(bar,'(V)')
axis square
xlabel('x (\mum)')
ylabel('y (\mum)');

Fig2=subplot(1,2,2);
imagesc(x,y,PhaseTauMask*180/pi)
title('E Phase')
colormap(Fig2,PhaseColormap)
bar=colorbar;
xlabel(bar,'(deg)')
axis square
caxis([-180 180])
xlabel('x (\mum)')
ylabel('y (\mum)')
saveas(F3,strcat(dir,'\','3.Tau.png'))

% 
% %% NEAR FIELD IN COMPLEX PLANE VISUALIZATION - NON ROTATED
% 
% [X1,Y1,X1rot,Y1rot,X1subs,Y1subs,X1subsrot,Y1subsrot,MaxEigenVect1,MinEigenVect1,SideB1xrotNF,SideB1yrotNF]=Covariance(K,C2*SideB2x,-C1*SideB1x);
% 
% FTau1=figure('units','normalized','outerposition',[0 0 1 1]);
% scatter(X1subs,Y1subs,'.','b')
% hold on
% scatter(X1,Y1,'.','r')
% quiver(0,0,MaxEigenVect1(1),MaxEigenVect1(2),'-m','LineWidth',2);
% quiver(0,0,MinEigenVect1(1),MinEigenVect1(2), '-g', 'LineWidth',2);
% hold off
% axis equal
% title('NEAR-FIELD (X COMPONENT - METHOD 1)')
% legend('substate','patch')
% xlabel('S_2_,_2 (V)')
% ylabel('S_2_,_1 (V)')
% saveas(FTau1,strcat(dir,'\','TauScatter.png'))
% %% NEAR FIELD IN COMPLEX PLANE VISUALIZATION - ROTATED
% 
% [X1,Y1,X1rot,Y1rot,X1subs,Y1subs,X1subsrot,Y1subsrot,MaxEigenVect1,MinEigenVect1,SideB1xrotNF,SideB1yrotNF]=Covariance(K,C2*SideB2xrot,-C1*SideB1xrot);
% 
% FTau2=figure('units','normalized','outerposition',[0 0 1 1]);
% scatter(X1subs,Y1subs,'.','b')
% hold on
% scatter(X1,Y1,'.','r')
% quiver(0,0,MaxEigenVect1(1),MaxEigenVect1(2),'-m','LineWidth',2);
% quiver(0,0,MinEigenVect1(1),MinEigenVect1(2), '-g', 'LineWidth',2);
% hold off
% axis equal
% title('NEAR-FIELD (X COMPONENT - METHOD 2)')
% legend('substate','patch')
% xlabel('S_2_,_2 (V)')
% ylabel('S_2_,_1 (V)')
% saveas(FTau2,strcat(dir,'\','TauScatterRotated.png'))

%% Histogram

% Phase

PhaseStruc=PhaseTau(StrucInd)*180/pi;
PhaseSubst=PhaseTau(SubstInd)*180/pi;

angle=[-180:360/36:180];

F4=figure('units','normalized','outerposition',[0 0 1 1]);
% suptitle('NEAR-FIELD (METHOD 1)')
Fig1=subplot(1,2,1);
histogram(PhaseStruc,angle);
hold on
histogram(PhaseSubst,angle);
title('Phase Histogram')
axis square
legend({'structure','substrate'})
xlabel('phase (deg)')
ylabel('Ocurrences');

Fig2=subplot(1,2,2);
imagesc(x,y,PhaseTau*180/pi)
title('E Phase')
colormap(Fig2,PhaseColormap)
bar=colorbar;
xlabel(bar,'(deg)')
axis square
caxis([-180 180])
xlabel('x (\mum)')
ylabel('y (\mum)')
saveas(F4,strcat(dir,'\','4.Phase.png'))

%% Phase Evolution

Nframes=101;
PhaseSpace=linspace(0,2*pi,Nframes);
dt=linspace(0,1,Nframes);

PhaseMovie(Nframes)=struct('cdata',[],'colormap',[]);

F5=figure();
for jj=1:Nframes
    PhaseEvol=PhaseTau+ones(N).*PhaseSpace(jj);
    Title=strcat('Phase (t=',num2str(dt(jj)),' period)');
    imagesc(x,y,wrapToPi(PhaseEvol.*MaskX)*180/pi,[-180, 180]);
    axis square; axis image; 
    colormap(PhaseColormap); colorbar
    xlabel('x (\mum)')
    ylabel('y (\mum)')
    title(Title)
%     set(gca,'FontSize',fst);
    PhaseMovie(jj)=getframe;
%     pause(0.01)
    
end

% Save video

v=VideoWriter('phase_evol_bowtie_4.avi');
v.FrameRate=10;
open(v);
writeVideo(v,PhaseMovie);
close(v);