close all

dt=0.01;

t=-1:dt:1;
M=size(t);
M=M(2);

u=linspace(-1/(2*dt),1/(2*dt),M);

f1=10;
f2=5;

s1=2*cos(2*pi*f1*t);
s2=2*cos(2*pi*f2*t);

s=s1.*s2;

figure(1)
plot(t,s1,t,s2,t,s)

figure(2)
plot(u,abs(fftshift(fft(s))/(M/2)))
