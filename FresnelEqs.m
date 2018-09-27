clear all
close all
clc

theta1=45*pi/180;
n1=1;
n2=2.4028;      % ZnSe at 10.6um

theta2=asin(n1*sin(theta1)/n2);

rp=(n1*cos(theta2)-n2*cos(theta1))/(n1*cos(theta2)+n2*cos(theta1));

rs=(n1*cos(theta1)-n2*cos(theta2))/(n1*cos(theta1)+n2*cos(theta2));

betacalc=atan2(rp,rs)*180/pi;

%% For a given beta

beta=30*pi/180;

Ei=1;   Eip=Ei*cos(beta);   Eis=Ei*sin(beta);

Erp=rp*Eip;     Ers=rs*Eis;

alpha=atan2(Ers,Erp)*180/pi

alphaprime=180-abs(alpha)
