function [Patch1Mean,Patch2Mean,Patch1Std,Patch2Std,Patch1CV,Patch2CV,MeanStd]=MeanAndStd(Patches1,Patches2)

global Dy Dx frame Xtop Xbottom Xleft Xright Ytop Ybottom Yleft Yright

% Select and store individual patches

SinglePatch1=zeros(Dx,Dy,4);
SinglePatch2=zeros(Dx,Dy,4);

for j=1:4
    SinglePatch1(:,:,j)=Patches1(Xtop(j):Xbottom(j),Xleft(j):Xright(j));
    SinglePatch2(:,:,j)=Patches2(Ytop(j):Ybottom(j),Yleft(j):Yright(j));
end

% Mean between patches

Patch1Mean=mean(SinglePatch1,3);
Patch2Mean=mean(SinglePatch2,3);

% Standar dev between patches

Patch1Std=std(SinglePatch1,0,3);
Patch2Std=std(SinglePatch2,0,3);

% Coefficient of Variation between patches (= std/mean)

Patch1CV=Patch1Std./abs(Patch1Mean);

N=size(Patch1CV);
N=N(2);

% for i=1:N
%     for j=1:N
%         if Patch1CV(i,j)>100
%             Patch1CV(i,j)=100;
%         end
%     end
% end

Patch2CV=Patch2Std./abs(Patch2Mean);

% for i=1:N
%     for j=1:N
%         if Patch2CV(i,j)>100
%             Patch2CV(i,j)=100;
%         end
%     end
% end

% Standard dev mean

MeanSideB1xStd=mean(mean(Patch1Std));
MeanSideB1yStd=mean(mean(Patch2Std));

MeanStd=[MeanSideB1xStd MeanSideB1yStd];
