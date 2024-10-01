clear all

N=6;J=-1;hz=-0;hx=-0.4;perBC=1;
CD={'k','b','r','g','m','c','y'};
pool={'I','Y','ZY','YZ'};
FalqonCD=2;% 0 - Original FALQON , 1 - FALOQN with Y, 2 - FALOQN with YZ
Colorcount=0;
beta=0;gamma=1;alpha=1;
HP=HamNN(N,3,3,J,perBC)+HamOnsite(N,3,hz)+HamOnsite(N,1,hx);
[EV,D]=eig(HP);

fenergy= sprintf('dataCD/MFIM/FalqonCDenergy%d.%f.%f.%f.%d.mat',N,J,hz,hx,FalqonCD);
fBeta= sprintf('dataCD/MFIM/FalqonCDBeta%d.%f.%f.%f.%d.mat',N,J,hz,hx,FalqonCD);
fGamma= sprintf('dataCD/MFIM/FalqonCDGamma%d.%f.%f.%f.%d.mat',N,J,hz,hx,FalqonCD);

E=importdata(fenergy);
B=importdata(fBeta);
G=importdata(fGamma);

Tlength=200;
energy=E(1:Tlength);
Beta=B(1:Tlength)
Gamma=G(1:Tlength)
dt=0.01;T=0:dt:dt*length(energy)-dt;


% Plots

subplot(1,3,1)
plot(T,real(energy(:,1))-D(1,1),'Color',CD{FalqonCD+Colorcount+1},'LineWidth',3)
hold on
ax1=gca;
xlabel('Time')
ylabel('Energy')
hold on
axis([0 length(T)*dt 0 real(energy(1,1))-D(1,1)])
fontsize(24,"points")

subplot(1,3,2)
plot(T,real(Beta(:,1)),'Color',CD{FalqonCD+Colorcount+1},'LineWidth',3)
hold on
%plot(T,real(betadot(:,1)),'-o','Color',CD{FalqonCD+1})
ax2=gca;
xlabel('Time')
ylabel('\beta(t)')
hold on
%axis([0 length(T)*dt min(real(Beta)) max(real(Beta))])
fontsize(24,"points")
str = sprintf('dt=%.2f,  J=%d, h_z=%.2f, h_x=%.2f',dt,J,hz,hx);
title(str)

subplot(1,3,3)
plot(T,real(Gamma(:,1)),'Color',CD{FalqonCD+Colorcount+1},'LineWidth',3)
hold on
%axis([0 length(T)*dt min(real(Gamma)) max(real(Gamma))])
fontsize(24,"points")
%plot(T,real(betadot(:,1)),'Color',CD{FalqonCD+2})
ax3=gca;
xlabel('Time')
ylabel('\gamma(t)')
