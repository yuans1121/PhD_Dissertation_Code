function [X,Y,Xrot,Yrot,Xsubs,Ysubs,Xsubsrot,Ysubsrot,MaxEigenVect,MinEigenVect,SideBxrot,SideByrot]=Covariance(K,SideBx,SideBy)

global Mx My MaskX MaskY

X=zeros(1,K);
Y=X;

Xsubs=zeros(1,Mx*My-K);
Ysubs=Xsubs;

k=1;
ksubs=1;
for ii=1:Mx
    for jj=1:My
        if MaskX(jj,ii)==1 && MaskY(jj,ii)==1
            X(k)=SideBx(jj,ii);
            Y(k)=SideBy(jj,ii);
            k=k+1;
        else
            Xsubs(ksubs)=SideBx(jj,ii);
            Ysubs(ksubs)=SideBy(jj,ii);
            ksubs=ksubs+1;
        end
    end
end

%% Center the distributions

% Xmean=mean(X);
% Ymean=mean(Y);
% 
% X=X-Xmean;
% Y=Y-Ymean;
% 
% Xsubsmean=mean(Xsubs);
% Ysubsmean=mean(Ysubs);
% 
% Xsubs=Xsubs-Xsubsmean;
% Ysubs=Ysubs-Ysubsmean;

%% Covariance Matrix

CovSB=cov(X,Y);
[EigenVect,EigenVal]=eig(CovSB);

% Find max eigenval & eigenvect

MaxEigenVal=max(max(EigenVal));
[MaxEigenVectInd,r]=find(EigenVal==MaxEigenVal);
MaxEigenVect=EigenVect(:,MaxEigenVectInd);

% Find min eigenval & eigenvect

if MaxEigenVectInd==1
    MinEigenVal=max(EigenVal(:,2));
    MinEigenVect=EigenVect(:,2);
else
    MinEigenVal=max(EigenVal(:,1));
    MinEigenVect=EigenVect(:,1);
end
    
% Angle between max eigenval & x-axis

angle=atan2(MaxEigenVect(2),MaxEigenVect(1));
angledeg=angle*180/pi

% Counterclockwise rotation

R=[cos(angle) sin(angle);-sin(angle) cos(angle)];

XYrot=R*[X;Y];
Xrot=XYrot(1,:);
Yrot=XYrot(2,:);
XYsubsrot=R*[Xsubs;Ysubs];
Xsubsrot=XYsubsrot(1,:);
Ysubsrot=XYsubsrot(2,:);

MaxEigenVect=MaxEigenVect*sqrt(MaxEigenVal);
MinEigenVect=MinEigenVect*sqrt(MinEigenVal);

% t=linspace(0,2*pi,50);
% Xellip=xSemiAxis*cos(t);
% Yellip=ySemiAxis*sin(t);

%% Sidebands rotation

SideBxrot=zeros(My,Mx);
SideByrot=SideBxrot;

k=1;
ksubs=1;
for ii=1:Mx
    for jj=1:My
        if MaskX(jj,ii)~=0 && MaskY(jj,ii)~=0
            SideBxrot(jj,ii)=Xrot(k);
            SideByrot(jj,ii)=Yrot(k);
            k=k+1;
        else
            SideBxrot(jj,ii)=Xsubsrot(ksubs);
            SideByrot(jj,ii)=Ysubsrot(ksubs);
            ksubs=ksubs+1;
        end
    end
end

end