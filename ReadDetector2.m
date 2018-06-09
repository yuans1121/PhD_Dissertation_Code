clear all
close all

Folder1='Data\2018-02-07\';

Folder2='Results\2018-02-07\';
mkdir(Folder2);

Signal='1';
WindowSize=2000; % um
z='755';
Series='4';
Step=50; % um
 
Size=num2str(WindowSize);

Scan=strcat('Detector_Signal_',Signal,'_',Size,'um_z_',z,'in_',Series);

f='.lvm';

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

Title=strcat('AFM Signal ',Signal,' - Detector at z=0.',z,' in [WindowSize=',num2str(WindowSize),'\mum, Step=',num2str(Step),'\mum]');
% Title=strcat('AFM + Reference Signal ',Signal,' - Detector at z=0.',z,' in [WindowSize=',num2str(WindowSize),'\mum, Step=',num2str(Step),'\mum]');

Fig=figure('units','normalized','outerposition',[0 0 1 1]);
% Fig=figure(1);
subplot(1,2,1)
suptitle(Title)
mesh(X,Y,abs(Z))
xlabel('x (\mum)')
ylabel('y (\mum)')
zlabel('Signal')
colormap(hot)
% colorbar
title('Detector Signal 3D')
axis square

% figure(2)
subplot(1,2,2)
pcolor(X,Y,abs(Z))
xlabel('x (\mum)')
ylabel('y (\mum)')
colormap(hot)
colorbar
title('Detector Signal 2D')
axis square

saveas(Fig,strcat(Folder2,'Signal_2ndHarmonic_',Signal,'_',Size,'um_',z,'_',Series,'.png'))