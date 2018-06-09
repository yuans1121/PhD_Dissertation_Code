%% ALIGN
% This program loads images taken from SNOM and finds the aplitude and 
% phase of the scattered near-field. 
% 1. Channels to be considered for analysis are specified in kk2h
% 2. By ginput, select a reference pixel for each image. All images will be
%    aligned using this reference
% 3. 

%% LOAD & SET REFERENCE ON THE IMAGES

clear all
close all
%load 2ndharmonic_dev2_18may17 
load 2016-07-26_squares_2and3H

% start=4;        % Channel start for 2nd harmonic FORWARD
start=3;        % Channel start for 2nd harmonic FORWARD

kk2h=start:3:3*12;      % Selects 2nd harmonic for 12 diff positions
kk3h=kk2h+1;                % Selects 3nd harmonic for 6 diff positions
kktopo=kk2h-1;              % Selects topography for 6 diff positions

mm2h=R(:,:,kk2h);      % 3D matrix contains the 6 images

deltax=x(2)-x(1);
deltay=y(2)-y(1);

nfiles=size(mm2h,3);        % Number of images to be analyse
nx=size(mm2h,1);                
ny=size(mm2h,2);

for kk=1:nfiles
    mm2h(:,:,kk)=flipud(mm2h(:,:,kk));       % Rotates matrix 90 degrees
end



for kk=1:nfiles;                        % Select a reference pixel via ginput using cursor
    figure(1);imagesc(mm2h(:,:,kk));axis square; axis xy
    [xref(kk),yref(kk)]=ginput;
end

%% IMAGES ALIGNMENT

fslabel=18;
fstitle=24;
fsticks=16;

iirefref=xref(1);      % Reference pixel
jjrefref=yref(1);      % Reference pixel

iirel=round(xref-iirefref);
jjrel=round(yref-jjrefref);

maxabs=max([max(abs(iirel)) max(abs(jjrel))]);

border=maxabs+10;

Nx=nx+2*border;
Ny=ny+2*border;

mm_align=zeros(Ny,Nx,nfiles);
mask=zeros(Ny,Nx,nfiles);
mm_align2=zeros(Ny,Nx,nfiles);

for kk=1:nfiles
    fromii=border-iirel(kk)+1;
    untilii=fromii+nx-1;
    fromjj=border-jjrel(kk)+1;
    untiljj=fromjj+ny-1;
    mm_align(fromjj:untiljj,fromii:untilii,kk)=mm2h(:,:,kk);
    mask(fromjj:untiljj,fromii:untilii,kk)=1;
    
    figure(5);imagesc(mm_align(:,:,1));axis square; axis xy;
    figure(6); hold off
    imagesc(mm_align(:,:,kk));axis square; axis xy;
    hold on;
%     iipuntorefkk=desdeii+iirefref+iirel(kk);
%     jjpuntorefkk=desdejj+jjrefref+jjrel(kk);
%     iipuntoref=desdeii+iirefref;
%     jjpuntoref=desdejj+jjrefref;
%     plot(iipuntorefkk,jjpuntorefkk,'ws','LineWidth',3);
%     hold on 
%     line([iipuntoref iipuntoref],[1 174],'Color','w');
%     line([1 174],[jjpuntoref jjpuntoref],'Color','w');
%     line([iipuntorefkk iipuntorefkk],[1 174],'Color','r');
%     line([1 174],[jjpuntorefkk jjpuntorefkk],'Color','r');

    
    %plot(iipuntoref,jjpuntoref,'wo','LineWidth',3);
    Text=strcat('iirel=',num2str(iirel(kk)),' jjrel=',num2str(jjrel(kk)));
    text(10,160,Text,'FontSize',fslabel,'Color','w');
    hold off      
    %subplot(2,2,3); imagesc(mascara(:,:,1));axis square; axis xy;
    %subplot(2,2,4); imagesc(mascara(:,:,kk));axis square; axis xy;
    pause(0.2)
end

maskprod=prod(mask,3);
[zerosii,zerosjj]=find(maskprod==0);

for kk=1:nfiles
    mm_align2(:,:,kk)=mm_align(:,:,kk).*maskprod;
    figure(2);imagesc(mm_align2(:,:,kk));axis square; axis xy;
    pause(0.2)
end

line1=maskprod(round(size(maskprod,1)/2),:);
dline1=diff(line1);        % Takes the "derivative"
INDEX1=find(dline1~=0);
index1x=INDEX1(1)+1;
index2x=INDEX1(2);
sizex=index2x-index1x+1;
line2=maskprod(:,round(size(maskprod,1)/2));
dline2=diff(line2);
INDEX2=find(dline2~=0);
index1y=INDEX2(1)+1;
index2y=INDEX2(2);
sizey=index2y-index1y+1;

mm_align3=zeros(sizey,sizex,nfiles);

xx=[0:deltax:deltax*(sizex-1)];
yy=[0:deltay:deltay*(sizey-1)];

for kk=1:nfiles
    mm_align3(:,:,kk)=mm_align2(index1y:index2y,index1x:index2x,kk);
    figure(3);imagesc(xx,yy,mm_align3(:,:,kk));axis equal; axis xy; axis image
    pause(0.5)
end

%% some analysis

iilineasustrato=110;
lineasustrato=mm_align3(iilineasustrato,:,:);
lineasustrato=squeeze(lineasustrato);

iilinea=20;
lambda0=10.2;
paso=0.5;

figure(3);
hold off
imagesc(mm_align3(:,:,1));axis equal; axis xy; axis image
hold on
line([1 size(mm_align3,2)],[iilinea iilinea] ,'Color','w')
hold off

linea=mm_align3(iilinea,:,:);
linea=squeeze(linea);
deltafase=360*2*paso/lambda0;
fases=[0:deltafase:deltafase*(nfiles-1)];

figure(4);
imagesc(fases,xx,sqrt(linea));

figure(5);
imagesc(fases,xx,sqrt(linea)./sqrt(lineasustrato));

%% analysis based on ed kinzel scripts.

posiciones=[25.00 25.50 26.00 26.50 27.00 27.50 28.00 28.50 29.00 29.50 30.00 30.50];
q=(posiciones-posiciones(1))./10.6;

[n1,n2,nfiles]=size(mm_align3);

posicion=[0:paso:(nfiles-1)*paso];
k0=4*pi/10.6;   % propagation constant
%q=posicion;
A=[ones(length(q),1) cos(q*k0)' sin(q*k0)'];                                % form least squares matrix

for kk=1:nfiles
z(:,:,kk) = mm_align3(:,:,kk)- mean(mean(mm_align3(:,:,kk)));
end


for u=1:n1,                                                                  % step over x
    for v=1:n2,                                                              % step over y
        b=permute(z(u,v,:),[3,1,2]);                                        % form signal vector as function of displacement
        ia=(A'*A)\A'*b;                                                     % project onto fitted space
        C4(u,v)=ia(1);                                                      % store constant value
        C5(u,v)=ia(2);                                                      % store cosine constant
        C6(u,v)=ia(3);                                                      % store sine constant
    end
end
aE=sqrt(C5.^2+C6.^2);                                                       % calculate amplitude
% aE=imrotate(aE,-1,'bilinear','crop')% Rotates image and crops to fit
% original size
pE=atan2(C6,C5);                                                            % calculate phase

% unwrapping

% for v=1:n2,                                                                  % step over lines (step over y)
%     for u=1:n1,
%             if sign(pE(u,v))==1,
%         pE(u,v)=rem(pE(u,v),2*pi);
%             elseif sign(pE(u,v))==-1
%    pE(u,v)=rem(pE(u,v),2*pi)+2*pi;
%             end
%     end
% end

%% figuras kinzel

mz=max(max(aE));                                                            % maximum amplitude

rm=(0:31)'/32; gm=[rm; 1; flipud(rm)]; rm=[rm; ones(33,1)]; bm=flipud(rm);  % bulid red-white-blue color map for phase
cm1=[rm gm bm];                                                             
gm=[zeros(32,1); (0:31)'/31]; rm=[(0:31)'/31; ones(32,1)]; bm=zeros(64,1);  % build black-red-yellow color map for amplitude
cm2=[rm gm bm];                                                             

f1=figure(11); set(f1,'name','Amplitude','NumberTitle','off')                   % plot signal amplitude
%srf2=surf(xx,yy,aE,'linestyle','none','facelighting','phong');
%srf2=surf(xx,yy,aE,'linestyle','none','facelighting','phong');
srf=imagesc(xx,yy,aE)
% surf(x,y,z) plots surface plot
%box on; grid off; view(2); 
axis equal; axis tight; axis xy
colormap(cm2); colorbar; 
%set(gca,'fontsize',18,'Ytick',[-6,-4,-2,0,2,4,6],'Xtick',[-6,-4,-2,0,2,4,6])
set(gca,'fontsize',18)%,'Ytick',[-6,-4,-2,0,2,4,6],'Xtick',[-6,-4,-2,0,2,4,6])
xlabel('\itx \rm(\mum)','fontweight','bold','FontName','Helvetica'); 
ylabel('\ity \rm(\mum)','fontweight','bold','FontName','Helvetica'); 
title('Amplitude','FontSize',fslabel+4)
colorbar('box','on')

%dlmwrite('phase.txt',pE', 'delimiter','\t');                                  % export .txt file of amplitude data in matrix and name ampl.txt (will put in same folder as matlab program)
%type phase.txt;

%dlmwrite('ampl.txt',aE, 'delimiter','\t');                                  % export .txt file of amplitude data in matrix and name ampl.txt (will put in same folder as matlab program)
%type ampl.txt;

%dlmwrite('mirrorpositions.txt',q, 'delimiter','\t');                       % export .txt file of reference mirror positions and name mirror positions.txt (will put in same folder as matlab program)
%type mirrorpositions.txt;

f2=figure(12); 
set(f2,'name','Phase','NumberTitle','off')                       % plot signal phase
%srf=surf(xx,yy,pE,'linestyle','none','facelighting','phong');
srf=imagesc(xx,yy,pE)
%box on; grid off; %view(2); 
axis equal; axis tight; axis xy 
colormap(cm1); colorbar;
xlabel('x [\mum]'); ylabel('y [\mum]');
set(gca,'fontsize',18)
title('Phase','FontSize',fslabel+4)


f3=figure(13); set(f3,'name','DC Residual','NumberTitle','off')                 % plot DC residual
%srf=surf(xx,yy,C4,'linestyle','none','facelighting','phong');
srf=imagesc(xx,yy,C4);
%box on; grid off; view(2); 
axis equal; axis tight; axis xy
colormap(cm1); colorbar;
xlabel('x [\mum]'); ylabel('y [\mum]');
set(gca,'fontsize',18)
title('DC Residual','FontSize',fslabel+4)
