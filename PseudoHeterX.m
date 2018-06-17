clear all
close all
clc

% Initialization

lambda=10.6; %microns

theta=60*pi/180;  % SNOM incidence angle

Vrms=0.106; % Volts

% MaxTau=50.4179; % Good diodes: From 2017-11-21_diode_h_5_PH_X_3
% MaxTau=26.2718; % Bad diodes diodes: From 2017-12-12_diode_g_2_PH_X_2
MaxTau=1;
% MaxTau=1;

run=1;

%% LOAD
    
load 2018-04-07_leaky_E_4_V_PH_X_1

dir1='Presentation';
dir2='2018-04-07';
dir3='Leaky_E4_V_1';
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
thrX=0.40*(MaxTopo-MinTopo)+MinTopo;

StrucInd=find(TopoX>thrX);
SubstInd=find(TopoX<thrX);

% [MaskX,Kx]=Mask(TopoX,thrX); 

MaskX=zeros(N,N);
MaskX(StrucInd)=1;
MaskX(SubstInd)=0;

% MaskX(1:12,:)=0;
% MaskX(:,1:16)=0;

Square=ones(3);
MaskX=conv2(MaskX,Square,'same');

StrucInd=find(MaskX>0);
SubstInd=find(MaskX==0);
MaskX(StrucInd)=1;

    F1=figure('units','normalized','outerposition',[0 0 1 1]);

    % X component

    subplot(1,2,1);
    imagesc(x,y,TopoX)
    % imagesc(x,y(1:N/2),TopoX(1:N/2,:))
    % imagesc(x,y,rot90(TopoX,3))
    title('Topography X')
    colormap(gray(256))
    % axis square
    xlabel('x (\mum)')
    ylabel('y (\mum)')
    daspect([1 1 1])

    subplot(1,2,2)
    imagesc(x,y,MaskX)
    % imagesc(x,y(1:N/2),MaskX(1:N/2,:))
    % imagesc(MaskX(1:N/2,:))
    % imagesc(x,y,rot90(MaskX,3))
    title('Mask X')
    colormap(gray(256))
    % axis square
    xlabel('x (\mum)')
    ylabel('y (\mum)')
    daspect([1 1 1])
    saveas(F1,strcat(dir,'\','1.Topography.png'))

if run==1

    %% ANALYSIS

    % 1ST SIDEBAND
    SideB1x=M(:,:,i+2);  % X component

    %2ND SIDEBAND
    SideB2x=M(:,:,i+3);  % X component

    % Phase color map (as used by Hillenbrand et al.)

    up=(1:89)'/90;
    dw=flipud(up);
    Zeros=zeros(89,1);
    Ones=ones(89,1);
    rm=[0 ; Zeros ; 0 ; up ; 1 ; Ones ; 1; dw];
    gm=[0 ; Zeros; 0 ; up ; 1 ; dw ; 0 ; Zeros]; 
    bm=[0 ; up ; 1 ; Ones ; 1 ; dw ; 0 ; Zeros]; 
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
        % imagesc(x,y(1:N/2),SideB1(1:N/2,:))
        % imagesc(x,y,rot90(SideB1,3))
        title('1st Sideband')
        colormap(FigSB1,hot)
        bar=colorbar;
        xlabel(bar,'(V)')
        axis square
        xlabel('x (\mum)')
        ylabel('y (\mum)')
        daspect([1 1 1])

        % 2ND sideband

        FigSB2=subplot(1,2,2);
        imagesc(x,y,SideB2)
        % imagesc(x,y(1:N/2),SideB2(1:N/2,:))
        % imagesc(x,y,rot90(SideB2,3))
        title('2nd Sideband')
        colormap(FigSB2,hot)
        bar=colorbar;
        xlabel(bar,'(V)')
        axis square
        xlabel('x (\mum)')
        ylabel('y (\mum)')
        daspect([1 1 1])
        saveas(F2,strcat(dir,'\','2.SideBands.png'))

    %% NEAR FIELD

    % Gamma calculation

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

    ModulTau=abs(Tau)./MaxTau;
    % MaxTau=max(max(ModulTau))
    % ModulTau(ModulTau>60)=0;
    PhaseTau=atan2(imag(Tau),real(Tau));

    Y=ones(My,Mx);
    for ii=1:1:Mx
       Y(:,ii)=y; 
    end
    PhaseCorr=PhaseTau-(2*pi/lambda)*Y*sin(theta);
    PhaseCorrW=wrapToPi(PhaseCorr);

    PhaseTauMask=PhaseTau.*MaskX;
    PhaseCorrWMask=PhaseCorrW.*MaskX;
    
        % Plot 1: Ez vs E Phase

        F31=figure('units','normalized','outerposition',[0 0 1 1]);
        Fig1=subplot(1,2,1);
        imagesc(x,y,ModulTau)
        title('|Ez|')
        colormap(Fig1,parula)
        bar=colorbar;
        % caxis([0 1])
        xlabel(bar,'(V)')
        axis square
        xlabel('x (\mum)')
        ylabel('y (\mum)');
        daspect([1 1 1])

        Fig2=subplot(1,2,2);
        imagesc(x,y,PhaseTauMask)
        title('Phase (Raw)')
        colormap(Fig2,PhaseColormap)
        bar=colorbar;
        xlabel(bar,'(rad)')
        axis square
        caxis([-pi pi])
        xlabel('x (\mum)')
        ylabel('y (\mum)')
        daspect([1 1 1])
        saveas(F31,strcat(dir,'\','3.Ez_vs_Phase.png'))

        % Plot 2: Ez vs E Phase (Corrected)

        F32=figure('units','normalized','outerposition',[0 0 1 1]);
        Fig1=subplot(1,2,1);
        imagesc(x,y,ModulTau)
        title('|Ez|')
        colormap(Fig1,parula)
        bar=colorbar;
        xlabel(bar,'(V)')
        axis square
        xlabel('x (\mum)')
        ylabel('y (\mum)');
        daspect([1 1 1])

        Fig2=subplot(1,2,2);
        imagesc(x,y,PhaseCorrWMask)
        title('Phase (Corrected)')
        colormap(Fig2,PhaseColormap)
        bar=colorbar;
        xlabel(bar,'(rad)')
        axis square
        caxis([-pi pi])
        xlabel('x (\mum)')
        ylabel('y (\mum)')
        daspect([1 1 1])
        saveas(F32,strcat(dir,'\','4.Ez_vs_PhaseCorr.png'))

        % Plot 3: Re{Ez} vs E Phase (Corrected)

        ReEz=ModulTau.*cos(PhaseCorrW);

        F33=figure('units','normalized','outerposition',[0 0 1 1]);
        Fig1=subplot(1,2,1);
        imagesc(x,y,ReEz)
        title('Re{E}')
        colormap(Fig1,parula)
        bar=colorbar;
        % caxis([0 1])
        xlabel(bar,'(V)')
        axis square
        xlabel('x (\mum)')
        ylabel('y (\mum)');
        daspect([1 1 1])

        Fig2=subplot(1,2,2);
        imagesc(x,y,PhaseCorrWMask)
        title('Phase (Corrected)')
        colormap(Fig2,PhaseColormap)
        bar=colorbar;
        xlabel(bar,'(rad)')
        axis square
        caxis([-pi pi])
        xlabel('x (\mum)')
        ylabel('y (\mum)')
        daspect([1 1 1])
        saveas(F33,strcat(dir,'\','5.ReEz_vs_PhaseCorr.png'))

        % Plot 4: Ez vs Re{E}

        F34=figure('units','normalized','outerposition',[0 0 1 1]);
        Fig1=subplot(1,2,1);
        imagesc(x,y,ModulTau)
        title('|Ez|')
        colormap(Fig1,parula)
        bar=colorbar;
        xlabel(bar,'(V)')
        axis square
        xlabel('x (\mum)')
        ylabel('y (\mum)');
        daspect([1 1 1])

        Fig2=subplot(1,2,2);
        imagesc(x,y,ReEz)
        title('Re{E}')
        colormap(Fig1,parula)
        bar=colorbar;
        xlabel(bar,'(V)')
        axis square
        xlabel('x (\mum)')
        ylabel('y (\mum)');
        daspect([1 1 1])
        saveas(F34,strcat(dir,'\','6.Ez_vs_ReEz.png'))

        % Plot 5: E Phase vs E Phase (Correcred)

        F35=figure('units','normalized','outerposition',[0 0 1 1]);
        Fig1=subplot(1,2,1);
        imagesc(x,y,PhaseTauMask)
        title('Phase (Raw)')
        colormap(Fig1,PhaseColormap)
        bar=colorbar;
        xlabel(bar,'(rad)')
        axis square
        caxis([-pi pi])
        xlabel('x (\mum)')
        ylabel('y (\mum)')
        daspect([1 1 1])

        Fig2=subplot(1,2,2);
        imagesc(x,y,PhaseCorrWMask)
        title('Phase (Corrected)')
        colormap(Fig2,PhaseColormap)
        bar=colorbar;
        xlabel(bar,'(rad)')
        axis square
        caxis([-pi pi])
        xlabel('x (\mum)')
        ylabel('y (\mum)')
        daspect([1 1 1])
        saveas(F35,strcat(dir,'\','7.Phase_vs_PhaseCorr.png'))

    %% NEAR FIELD IN COMPLEX PLANE VISUALIZATION - NON ROTATED
    
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

%     PhaseStruc=PhaseCorr(StrucInd);
%     PhaseSubst=PhaseCorr(SubstInd);
% 
%     angle=[-pi:2*pi/36:pi];
% 
%         F4=figure('units','normalized','outerposition',[0 0 1 1]);
%         % suptitle('NEAR-FIELD (METHOD 1)')
%         Fig1=subplot(1,2,1);
%         histogram(PhaseStruc,angle);
%         hold on
%         histogram(PhaseSubst,angle);
%         title('Phase Histogram')
%         axis square
%         legend({'structure','substrate'})
%         xlabel('phase (rad)')
%         ylabel('Ocurrences');
% 
%         Fig2=subplot(1,2,2);
%         imagesc(x,y,PhaseCorr)
%         % imagesc(x,y(1:N/2),PhaseTau(1:N/2,:))
%         % imagesc(x,y,rot90(PhaseTau,3))
%         title('E Phase')
%         colormap(Fig2,PhaseColormap)
%         bar=colorbar;
%         xlabel(bar,'(rad)')
%         axis square
%         caxis([-pi pi])
%         xlabel('x (\mum)')
%         ylabel('y (\mum)')
%         daspect([1 1 1])
%         % saveas(F4,strcat(dir,'\','4.Phase.png'))

    %% Phase Evolution

    Nframes=101;
    PhaseSpace=linspace(0,2*pi,Nframes);
    dt=linspace(0,1,Nframes);

    PhaseMovie(Nframes)=struct('cdata',[],'colormap',[]);

        F5=figure();
        for jj=1:Nframes
            PhaseEvol=PhaseCorrW+ones(N).*PhaseSpace(jj);
            Title=strcat('Phase (t=',num2str(dt(jj)),' period)');
            WrapPhase=wrapToPi(PhaseEvol.*MaskX);
            imagesc(x,y,WrapPhase,[-pi, pi]);
        %     imagesc(x,y(1:N/2),WrapPhase(1:N/2,:),[-pi, pi]);
        %     imagesc(x,y,rot90(WrapPhase,3),[-pi, pi]);
            axis square; axis image; 
            colormap(PhaseColormap); colorbar
            xlabel('x (\mum)')
            ylabel('y (\mum)')
        %     daspect([1 1 1])
            title(Title)
        %     set(gca,'FontSize',fst);
            drawnow
            PhaseMovie(jj)=getframe;
        %     pause(0.01)

        end

        % Save video

        v=VideoWriter(strcat(dir,'\','8.Phase_Evol.avi'));
        v.FrameRate=10;
        open(v);
        writeVideo(v,PhaseMovie);
        close(v);

    %% Histogram

    % Phase

    % ModulTauStruct=ModulTau(StrucInd);
    % SB2Struct=SideB2(SubstInd);
    % 
        % F6=figure('units','normalized','outerposition',[0 0 1 1]);
        % % suptitle('NEAR-FIELD (METHOD 1)')
        % Fig61=subplot(1,2,1);
        % histogram(ModulTau);
        % title('Sidenad 1 Histogram')
        % axis square
        % xlabel('Signal')
        % ylabel('Ocurrences');
        % 
        % Fig62=subplot(1,2,2);
        % histogram(SB2Struct);
        % title('Sidenad 2 Histogram')
        % axis square
        % xlabel('Signal')
        % ylabel('Ocurrences');
end