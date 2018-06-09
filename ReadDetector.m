clear all
close all


% Folder1='Data\2016-10-28\';
Folder1='Results\2017-02-02\';

Signal=1;           z=722;
Signal=num2str(Signal);  z=num2str(z);

% Scan=strcat('Detector_Signal_1000x1000um_',Signal,'_z_',z,'in');
% Scan=strcat('2ndHarmonic_Signal_300x300um_',Signal,'_z_',z,'in');
Scan=strcat('2ndHarmonic_Signal_300x300um_',Signal,'_z_',z,'in');


f='.lvm';

WindowSize=400; % um
Step=5; % um

N=floor(WindowSize/Step);

File=strcat(Folder1,Scan,f);
FileID=fopen(File);
data=fscanf(FileID,'%f\t%f\t%f\t%f\t%f',[5 inf]);
fclose(FileID);

x1=data(3,1); x2=data(3,N); X=x1:Step:x2;
y1=data(4,1); y2=data(4,N*N); Y=y1:Step:y2;

S=data(5,:);
Z=zeros(N,N);
for i=1:N
    for j=1:N
        Z(i,j)=S(N*(i-1)+j);
    end
end

Title=strcat('AFM + Reference Signal ',Signal,' - Detector at z=0.',z,' in [WindowSize=',num2str(WindowSize),'\mum, Step=',num2str(Step),'\mum]');

Fig=figure('units','normalized','outerposition',[0 0 1 1]);
% Fig=figure(1);
subplot(1,2,1)
suptitle(Title)
mesh(X,Y,Z)
xlabel('x (\mum)')
ylabel('y (\mum)')
zlabel('Signal')
colormap(hot)
% colorbar
title('Detector Signal 3D')
axis square

% figure(2)
subplot(1,2,2)
pcolor(X,Y,Z)
xlabel('x (\mum)')
ylabel('y (\mum)')
colormap(hot)
colorbar
title('Detector Signal 2D')
axis square

saveas(Fig,strcat(Folder1,'Signal_2ndHarmonic_',Signal,'_',z,'.png'))