clear all
close all
clc

% Initialization

lambda=10.6;    %microns

theta=60*pi/180;    % SNOM incidence angle

Vrms=0.104;     % Volts

run=1;  % run=0 stops after topography
        % run=1 continues after topography

savefigs=1; % savefigs=0 do not save images
            % savefigs=1 save images
            
savevid=0;  % savefigs=0 do not save video
            % savefigs=1 save video
        
saveMATLABfig=0;    % saveMATLABfig=0 do not save figure in matlab format
                    % saveMATLABfig=0 save figure in matlab format

%% LOAD
% Section loads the data and creates the respective folders

load 2018-10-22_sample_g_diode_11_PH_1

dir1='Analysis';
dir2='2018-10-22';
dir3='sample_g_diode_11_1';
dir4=strcat(dir2,'\',dir3);
mkdir(dir1,dir4);
dir=strcat(dir1,'\',dir4);

global Mx My MaskX MaskY
Mx=size(x);     My=size(y);
Mx=Mx(2);       My=My(2);

i=1; 


%% TOPOGRAPHY 
% Section creates the matrices for the topography maps for the X,Y
% components and its respective masks.

% X component

TopoX=M(:,:,i);
N=size(TopoX);
N=N(1);

MinTopo=min(min(TopoX));
MaxTopo=max(max(TopoX));
thrX=0.50*(MaxTopo-MinTopo)+MinTopo;

StrucInd=find(TopoX>thrX);
Kx=size(StrucInd); Kx=Kx(1);
SubstInd=find(TopoX<thrX); 

MaskX=zeros(N,N);
MaskX(StrucInd)=1;
MaskX(SubstInd)=0;

Square=ones(3);
MaskX=conv2(MaskX,Square,'same');

StrucInd=find(MaskX>0);
SubstInd=find(MaskX==0);
MaskX(StrucInd)=1;

% Y component 

TopoY=M(:,:,i+4);
N=N(1);

MinTopo=min(min(TopoY));
MaxTopo=max(max(TopoY));
thrY=thrX;

StrucInd=find(TopoY>thrY);
Ky=size(StrucInd); Ky=Ky(1);
SubstInd=find(TopoY<thrY); 

MaskY=zeros(N,N);
MaskY(StrucInd)=1;
MaskY(SubstInd)=0;

Square=ones(3);
MaskY=conv2(MaskY,Square,'same');

StrucInd=find(MaskY>0);
SubstInd=find(MaskY==0);
MaskY(StrucInd)=1;

if Kx<Ky
    K=Kx;
else
    K=Ky;
end

MaskX=MaskY;

    F1=figure(); % Topographies for the X & Y scans and its masks.

    % X component

    subplot(2,2,1)
    imagesc(x,y,TopoX)
    title('Topography X')
    colormap(gray(256))
%     axis square
    xlabel('x (\mum)')
    ylabel('y (\mum)')
    daspect([1 1 1])

    subplot(2,2,2)
    imagesc(x,y,MaskX)
    title('Mask X')
    colormap(gray(256))
%     axis square
    xlabel('x (\mum)')
    ylabel('y (\mum)');
    daspect([1 1 1])

    % Y component 

    subplot(2,2,3)
    imagesc(x,y,TopoY)
    title('Topography Y')
    colormap(gray(256))
%     axis square
    xlabel('x (\mum)')
    ylabel('y (\mum)');
    daspect([1 1 1])

    subplot(2,2,4)
    imagesc(x,y,MaskY)
    title('Mask Y')
    colormap(gray(256))
%     axis square
    xlabel('x (\mum)')
    ylabel('y (\mum)');
    daspect([1 1 1])
    pause(0.00001);
    frame_h = get(handle(gcf),'JavaFrame');
    set(frame_h,'Maximized',1);


if run==1

    %% ANALYSIS

    % 1ST SIDEBAND

    SideB1x=M(:,:,i+2);  % X component
    SideB1y=M(:,:,i+6);  % Y component

    %2ND SIDEBAND

    SideB2x=M(:,:,i+3);  % X component
    SideB2y=M(:,:,i+7);  % Y component

    % Phase color map (as used by Hillenbrand et al.)

    up=(1:89)'/90;
    dw=flipud(up);
    Zeros=zeros(89,1);
    Ones=ones(89,1);
    rm=[0 ; Zeros ; 0 ; up ; 1 ; Ones ; 1; dw];
    gm=[0 ; Zeros; 0 ; up ; 1 ; dw ; 0 ; Zeros]; 
    bm=[0 ; up ; 1 ; Ones ; 1 ; dw ; 0 ; Zeros]; 
    PhaseColormap=[rm gm bm];    

    %% SIDEBANDS - NON ROTATED

    % 1ST sideband
    SideB1=sqrt(SideB1x.^2+SideB1y.^2);
    SideB1phase=MaskX.*atan2(SideB1y,SideB1x)*180/pi;

    % 2ND sideband
    SideB2=sqrt(SideB2x.^2+SideB2y.^2);
    SideB2phase=MaskX.*atan2(SideB2y,SideB2x)*180/pi;

    
        % Display the X and Y components for the SBs, and the calculated modulus
        F2=figure();  

        % 1ST sideband
        % suptitle(strcat(dir3,' (Non Rotated)'))
        % suptitle('SIDEBANDS (METHOD 1)')
        FigSB1x=subplot(2,3,1);
        imagesc(x,y,SideB1x)
        title('1st Sideband X')
        colormap(FigSB1x,hot)
        bar=colorbar;
        xlabel(bar,'(V)')
%         axis square
        xlabel('x (\mum)')
        ylabel('y (\mum)')
        daspect([1 1 1])

        FigSB1y=subplot(2,3,2);
        imagesc(x,y,SideB1y)
        title('1st Sideband Y')
        colormap(FigSB1y,hot)
        bar=colorbar;
        xlabel(bar,'(V)')
%         axis square
        xlabel('x (\mum)')
        ylabel('y (\mum)')
        daspect([1 1 1])


        FigSB1=subplot(2,3,3);
        imagesc(x,y,SideB1)
        title('1st Sideband Modulus')
        colormap(FigSB1,hot)
        bar=colorbar;
        xlabel(bar,'(V)')
%         axis square
        xlabel('x (\mum)')
        ylabel('y (\mum)')
        daspect([1 1 1])

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
        bar=colorbar;
        xlabel(bar,'(V)')
%         axis square
        xlabel('x (\mum)')
        ylabel('y (\mum)')
        daspect([1 1 1])

        FigSB2y=subplot(2,3,5);
        imagesc(x,y,SideB2y)
        title('2nd Sideband Y')
        colormap(FigSB2y,hot)
        bar=colorbar;
        xlabel(bar,'(V)')
%         axis square
        xlabel('x (\mum)')
        ylabel('y (\mum)')
        daspect([1 1 1])

        FigSB2=subplot(2,3,6);
        imagesc(x,y,SideB2)
        title('2nd Sideband Modulus')
        colormap(FigSB2,hot)
        bar=colorbar;
        xlabel(bar,'(V)')
%         axis square
        xlabel('x (\mum)')
        ylabel('y (\mum)')
        daspect([1 1 1])
        
        % FigSB2ph=subplot(2,4,8);
        % imagesc(x,y,SideB2phase)
        % title('1st Sideband Phase')
        % colormap(FigSB2ph,rwb)
        % colorbar
        % axis square
        % xlabel('x (\mum)')
        % ylabel('y (\mum)')
        pause(0.00001);
        frame_h = get(handle(gcf),'JavaFrame');
        set(frame_h,'Maximized',1);

    %% NEAR FIELD (Sideband Modulus)

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
    
    % E-field
    
    Tau=C2*SideB2-1i*C1*SideB1;

    ModulTau=abs(Tau);
    PhaseTau=atan2(-imag(Tau),real(Tau));

    PhaseTauMask=PhaseTau.*MaskX;   % Raw Phase
    PhaseCorrW=PhaseCorW(y,PhaseTau,lambda,theta);
    PhaseCorrWMask=PhaseCorrW.*MaskX; % Corrected phase
    
        % Near-field magnitude and phase, using the SBs modulus 
        F3=figure();
        Fig1=subplot(1,2,1);
        imagesc(x,y,ModulTau)
        title('|Ez|')
        colormap(Fig1,parula)
        bar=colorbar;
        xlabel(bar,'(V)')
%         axis square
        xlabel('x (\mum)')
        ylabel('y (\mum)');
        daspect([1 1 1])

        Fig2=subplot(1,2,2);
        imagesc(x,y,PhaseCorrWMask)
        title('E Phase')
        colormap(Fig2,PhaseColormap)
        bar=colorbar;
        xlabel(bar,'(deg)')
%         axis square
        xlabel('x (\mum)')
        ylabel('y (\mum)')
        daspect([1 1 1])
        pause(0.00001);
        frame_h = get(handle(gcf),'JavaFrame');
        set(frame_h,'Maximized',1);

    %% NEAR FIELD (Only X component)

    Tau=C2*SideB2x-1i*C1*SideB1x;

    ModulTau=abs(Tau);
    PhaseTau=atan2(-imag(Tau),real(Tau));
    PhaseTauMask=PhaseTau.*MaskX;   % Raw Phase
    PhaseCorrW=PhaseCorW(y,PhaseTau,lambda,theta);
    PhaseCorrWMask=PhaseCorrW.*MaskX; % Corrected phase
    
        % Near-field magnitude and phase, using the SBs X component only
        F4=figure();
        Fig1=subplot(1,2,1);
        imagesc(x,y,ModulTau)
        title('|Ez|')
        colormap(Fig1,parula)
        bar=colorbar;
        xlabel(bar,'(V)')
%         axis square
        xlabel('x (\mum)')
        ylabel('y (\mum)');
        daspect([1 1 1])

        Fig2=subplot(1,2,2);
        imagesc(x,y,PhaseCorrWMask)
        title('E Phase')
        colormap(Fig2,PhaseColormap)
        bar=colorbar;
        xlabel(bar,'(deg)')
%         axis square
        xlabel('x (\mum)')
        ylabel('y (\mum)')
        daspect([1 1 1])
        pause(0.00001);
        frame_h = get(handle(gcf),'JavaFrame');
        set(frame_h,'Maximized',1);
        

    %% ROTATION

    [X1,Y1,X1rot,Y1rot,X1subs,Y1subs,X1subsrot,Y1subsrot,MaxEigenVect1,MinEigenVect1,SideB1xrot,SideB1yrot]=Covariance(K,SideB1x,SideB1y);

    [X2,Y2,X2rot,Y2rot,X2subs,Y2subs,X2subsrot,Y2subsrot,MaxEigenVect2,MinEigenVect2,SideB2xrot,SideB2yrot]=Covariance(K,SideB2x,SideB2y);

        % Scatter plots for the original SBs and rotated SBs
        F5=figure();
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
        daspect([1 1 1])

        subplot(2,2,2)
        scatter(X1rot,Y1rot,'.','r')
        % hold on
        % plot(Xellip1,Yellip1)
        % hold off
        axis equal
        title('1st Sideband (rotated)')
        xlabel('X (V)')
        ylabel('Y (V)')
        daspect([1 1 1])

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
        daspect([1 1 1])

        subplot(2,2,4)
        scatter(X2rot,Y2rot,'.','r')
        % hold on
        % plot(Xellip2,Yellip2)
        % hold off
        axis equal
        title('2nd Sideband (rotated)')
        xlabel('X (V)')
        ylabel('Y (Y)')
        daspect([1 1 1])
        pause(0.00001);
        frame_h = get(handle(gcf),'JavaFrame');
        set(frame_h,'Maximized',1);
        
    %% SIDEBANDS - ROTATED

    % 1ST sideband
    SideB1rot=sqrt(SideB1xrot.^2+SideB1yrot.^2);
    SideB1rotphase=MaskX.*atan2(SideB1yrot,SideB1xrot)*180/pi;

    % 2ND sideband
    SideB2rot=sqrt(SideB2xrot.^2+SideB2yrot.^2);
    SideB2rotphase=MaskX.*atan2(SideB2yrot,SideB2xrot)*180/pi;

    % 1ST sideband

        % Display the X and Y components for the rotated SBs, and the calculated modulus 
        F6=figure();
        % suptitle('SIDEBANDS (METHOD 2)')
        FigSB1xrot=subplot(2,3,1);
        imagesc(x,y,SideB1xrot)
        title('1st Sideband X')
        colormap(FigSB1xrot,hot)
        colorbar
        axis square
        bar=colorbar;
        xlabel(bar,'(V)')
        xlabel('x (\mum)')
        ylabel('y (\mum)');
        daspect([1 1 1])

        FigSB1yrot=subplot(2,3,2);
        imagesc(x,y,SideB1yrot)
        title('1st Sideband Y')
        colormap(FigSB1yrot,hot)
        bar=colorbar;
        xlabel(bar,'(V)')
%         axis square
        xlabel('x (\mum)')
        ylabel('y (\mum)');
        daspect([1 1 1])

        FigSB1rot=subplot(2,3,3);
        imagesc(x,y,SideB1rot)
        title('1st Sideband Modulus')
        colormap(FigSB1rot,hot)
        bar=colorbar;
        xlabel(bar,'(V)')
%         axis square
        xlabel('x (\mum)')
        ylabel('y (\mum)')
        daspect([1 1 1])

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
        bar=colorbar;
        xlabel(bar,'(V)')
%         axis square
        xlabel('x (\mum)')
        ylabel('y (\mum)')
        daspect([1 1 1])

        FigSB2yrot=subplot(2,3,5);
        imagesc(x,y,SideB2yrot)
        title('2nd Sideband Y')
        colormap(FigSB2yrot,hot)
        bar=colorbar;
        xlabel(bar,'(V)')
%         axis square
        xlabel('x (\mum)')
        ylabel('y (\mum)')
        daspect([1 1 1])

        FigSB2rot=subplot(2,3,6);
        imagesc(x,y,SideB2rot)
        title('2nd Sideband Modulus')
        colormap(FigSB2rot,hot)
        bar=colorbar;
        xlabel(bar,'(V)')
%         axis square
        xlabel('x (\mum)')
        ylabel('y (\mum)')
        daspect([1 1 1])

        % FigSB2rotph=subplot(2,4,8);
        % imagesc(x,y,SideB2rotphase)
        % title('1st Sideband Phase')
        % colormap(FigSB2rotph,rwb)
        % colorbar
        % axis square
        % xlabel('x (\mum)')
        % ylabel('y (\mum)')

        pause(0.00001);
        frame_h = get(handle(gcf),'JavaFrame');
        set(frame_h,'Maximized',1);

    %% NEAR FIELD - ROTATED

    Taurot=C2*SideB2xrot-1i*C1*SideB1xrot;

    ModulTaurot=abs(Taurot);
    PhaseTaurot=atan2(-imag(Taurot),real(Taurot));
    PhaseTaurotMask=PhaseTaurot.*MaskX;   % Raw Phase
    PhaserotCorrW=PhaseCorW(y,PhaseTaurot,lambda,theta);
    PhaserotCorrWMask=PhaserotCorrW.*MaskX; % Corrected phase
    
        % Near-field magnitude and phase, using the rotated SBs X component only
        F7=figure();
        % suptitle('NEAR-FIELD METHOD 2')
        Fig1=subplot(1,2,1);
        imagesc(x,y,ModulTaurot)
        title('|Ez| (Rotated)')
        colormap(Fig1,parula)
        bar=colorbar;
        xlabel(bar,'(V)')
%         axis square
        xlabel('x (\mum)')
        ylabel('y (\mum)');
        daspect([1 1 1])

        Fig2=subplot(1,2,2);
        imagesc(x,y,PhaserotCorrWMask)
        title('E Phase (Rotated)')
        colormap(Fig2,PhaseColormap)
        bar=colorbar;
        xlabel(bar,'(deg)')
        axis square
        xlabel('x (\mum)')
        ylabel('y (\mum)')
        daspect([1 1 1])
        pause(0.00001);
        frame_h = get(handle(gcf),'JavaFrame');
        set(frame_h,'Maximized',1);

    %% NEAR FIELD IN COMPLEX PLANE VISUALIZATION - NON ROTATED

    [X1,Y1,X1rot,Y1rot,X1subs,Y1subs,X1subsrot,Y1subsrot,MaxEigenVect1,MinEigenVect1,SideB1xrotNF,SideB1yrotNF]=Covariance(K,C2*SideB2x,-C1*SideB1x);

        % Near-field magnitude in the complex plane as a scatter plot.
        % Calculation was done using the original SBs x components
        FTau1=figure();
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
        pause(0.00001);
        frame_h = get(handle(gcf),'JavaFrame');
        set(frame_h,'Maximized',1);
        
    %% NEAR FIELD IN COMPLEX PLANE VISUALIZATION - ROTATED

    [X1,Y1,X1rot,Y1rot,X1subs,Y1subs,X1subsrot,Y1subsrot,MaxEigenVect1,MinEigenVect1,SideB1xrotNF,SideB1yrotNF]=Covariance(K,C2*SideB2xrot,-C1*SideB1xrot);

        % Near-field magnitude in the complex plane as a scatter plot.
        % Calculation was done using the rotated SBs x components
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
        pause(0.00001);
        frame_h = get(handle(gcf),'JavaFrame');
        set(frame_h,'Maximized',1);
        
    %% Save Figures
    if savefigs==1
        saveas(F1,strcat(dir,'\','1.Topography.png'))
        saveas(F2,strcat(dir,'\','2.SideBands.png'))
        saveas(F3,strcat(dir,'\','3.Ez_and_Phase-using_modulus.png'))
        saveas(F4,strcat(dir,'\','4.Ez_and_Phase-using_xcomp.png'))
        saveas(F5,strcat(dir,'\','5.SideBands_scatterplot.png'))
        saveas(F6,strcat(dir,'\','6.RotatedSideBands.png'))
        saveas(F7,strcat(dir,'\','7.Ez_and_Phase-using_rotSBs.png'))
    end
end