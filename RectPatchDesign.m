% Rectangular patch Design
% following Balanis pg 791

er=4.7447;      % relat permit. for ZnC
fr=2.83e13;     % (Hz) frequency at 10.6um
h=317e-9;       % (m) substrate thickness

c=3e8;          % (m/s) speed of light


W=(c/(2*fr))*sqrt(2/(er+1))

ereff=0.5*(er+1)+0.5*(er-1)/sqrt(1+12*h/W)

DeltaL=h*0.412*(ereff+0.3)*(W/h+0.264)/((ereff-0.258)*(W/h+0.8))

L=c/(2*fr*sqrt(ereff))-2*DeltaL