%Gamma calculation

lambda=10.6; % microns

Vrms=0.178; % Volts

V=sqrt(2)*Vrms;
scale=0.2218;
DeltaL=V/scale;
gamma=4*pi*DeltaL/lambda;

J1=besselj(1,gamma);
J2=besselj(2,gamma);

C1=1/J1;
C2=1/J2;