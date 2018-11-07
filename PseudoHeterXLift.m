clear all
close all
clc

% Initialization

lambda=10.6;    %microns

theta=60*pi/180;    % SNOM incidence angle

Vrms=0.104;     % Volts

%% LOAD
    
load 2018-10-01_lift_leaky_E4_PH_VLP_3

dir1='Analysis';
dir2='2018-10-01';
dir3='lift_leaky_E4_VLP_3';
dir4=strcat(dir2,'\',dir3);
mkdir(dir1,dir4);
dir=strcat(dir1,'\',dir4);

%% TOPOGRAPHY

y=MSB1(:,2,1);
n=size(y);  n=n(1);
nh=size(h); nh=nh(1);

TopoFwrd=flipud(MSB1(:,5,1));
TopoBwrd=flipud(MSB1(:,6,1));

z=TopoFwrd-min(TopoFwrd);

p=polyfit(y,z,1);
zfit=polyval(p,y);

zp=z-zfit+min(z);

% Y=zeros(1,6);
% Z=zeros(1,6);

% j=1;
% for i=3:1:n-3
%     if(TopoFwrd(i)<0.25)
%         if(TopoFwrd(i)<TopoFwrd(i-2) && TopoFwrd(i)<TopoFwrd(i+2))
%             Y(j)=y(i);
%             Z(j)=TopoFwrd(i);
%             j=j+1;
%         end
%     end
% end

% Y=Y-Y(1);
% Z=Z-Z(1);
% theta=atan(Z./Y);
% thetadeg=theta*180/pi;
% 
% theta=mean(theta(2:6));
% 
% TopoFwrdcorr=TopoFwrd-y*tan(theta);
% TopoBwrdcorr=TopoBwrd-y*tan(theta);
% 
% thr=0.20;
% 
% StrucInd=find(TopoFwrdcorr>thr);
% SubstInd=find(TopoFwrdcorr<thr);
% 
% Mask=zeros(1,n);
% Mask(StrucInd)=1;
% Mask(SubstInd)=0;

% diff=mean(abs(TopoFwrdcorr-TopoBwrdcorr));

angle=atan((z(1)-z(n))/y(n))*180/pi

%% ANALYSIS

m=size(h);  m=m(2);
% m=1;
SideB1=zeros(m,n);
SideB2=zeros(m,n);

for i=1:m
    SideB1(i,:)=MSB1(:,4,i)';
    SideB2(i,:)=MSB2(:,4,i)';
end

SideB1=flipud(SideB1);
SideB2=flipud(SideB2);

% Gamma calculation

V=sqrt(2)*Vrms;
scale=0.2218;
DeltaL=V/scale;
gamma=4*pi*DeltaL/lambda;

J1=besselj(1,gamma);
J2=besselj(2,gamma);

C1=1/J1;
C2=1/J2;

Vrms=(2.63*lambda/(4*pi))*scale/(sqrt(2));

Tau=C2*SideB2-1i*C1*SideB1;
    
MaxTau=max(abs(Tau));
ModulTau=abs(Tau);
ModulTauNorm=abs(Tau)/MaxTau;
    
PhaseTau=atan2(-imag(Tau),real(Tau));
% PhaseTauMask=PhaseTau.*Mask;
PhaseTauMask=PhaseTau;

Y=ones(nh,n);
for ii=1:1:nh
    Y(ii,:)=y; 
end

PhaseCorr=PhaseTau-(2*pi/lambda)*Y*sin(theta);
PhaseCorrW=wrapToPi(PhaseCorr);

% Phase color map (as used by Hillenbrand et al.)

up=(1:89)'/90;
dw=flipud(up);
Zeros=zeros(89,1);
Ones=ones(89,1);
rm=[0 ; Zeros ; 0 ; up ; 1 ; Ones ; 1; dw];
gm=[0 ; Zeros; 0 ; up ; 1 ; dw ; 0 ; Zeros]; 
bm=[0 ; up ; 1 ; Ones ; 1 ; dw ; 0 ; Zeros]; 
PhaseColormap=[rm gm bm];  

%% PLOTS

% h=flip(h);

% ModulTau=flip(ModulTau);
% PhaseTauMask=flip(PhaseTauMask);

Fig1=figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(1,2,1);
    plot(y,TopoFwrd,y,TopoBwrd,y,zfit)
    title('Topography')
    axis square
    xlabel('y (\mum)')
    ylabel('z (\mum)');

    subplot(1,2,2);
    plot(y,z,y,zp)
    title('Mask')
    axis square
    xlabel('y (\mum)')
    ylabel('z (\mum)');
        
Fig2=figure('units','normalized','outerposition',[0 0 1 1]);
    Fig21=subplot(2,1,1);
    imagesc(y,h,flipud(ModulTau))
    title('Field Modulus |Ez|')
    colormap(Fig1,parula)
    bar=colorbar;
    xlabel(bar,'(V)')
    xlabel('y (\mum)')
    ylabel('h (nm)');
    set(gca,'Ydir','normal')
    
    Fig22=subplot(2,1,2);
    imagesc(y,h,flipud(PhaseCorrW))
    title('Phase')
    colormap(Fig22,PhaseColormap)
    bar=colorbar;
    xlabel(bar,'(rad)')
    caxis([-pi pi])
    xlabel('y (\mum)')
    ylabel('h (nm)')
    set(gca,'Ydir','normal')
    pause(0.00001);
    frame_h = get(handle(gcf),'JavaFrame');
    set(frame_h,'Maximized',1);
    
Fig3=figure('units','normalized','outerposition',[0 0 1 1]);
    Fig21=subplot(2,1,1);
    imagesc(y,h,flipud(ModulTau))
    title('Field Modulus |Ez|')
    colormap(Fig1,parula)
    xlabel('y (\mum)')
    ylabel('h (nm)');
    set(gca,'Ydir','normal')
    
    Fig22=subplot(2,1,2);
    plot(y,z)
    title('Topography')
    xlim([0 15])
    xlabel('y (\mum)')
    ylabel('z (\mum)');
    pause(0.00001);
    frame_h = get(handle(gcf),'JavaFrame');
    set(frame_h,'Maximized',1);
    
Fig4=figure('units','normalized','outerposition',[0 0 1 1]);
    Fig21=subplot(2,1,1);
    imagesc(y,h,flipud(PhaseCorrW))
    title('Phase')
    colormap(Fig21,PhaseColormap)
    caxis([-pi pi])
    xlabel('y (\mum)')
    ylabel('h (nm)')
    set(gca,'Ydir','normal')
    
    Fig22=subplot(2,1,2);
    plot(y,z)
    title('Topography')
    xlim([0 15])
    xlabel('y (\mum)')
    ylabel('z (\mum)');
    pause(0.00001);
    frame_h = get(handle(gcf),'JavaFrame');
    set(frame_h,'Maximized',1);
    
% saveas(Fig2,strcat(dir,'\','1.modulus_phase.png'));
% saveas(Fig3,strcat(dir,'\','2.modulus.png'));
% saveas(Fig4,strcat(dir,'\','3.phase.png'));