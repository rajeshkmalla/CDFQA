clear all
%close all
tic
N=8;
J=1;
hz=0.809;
hx=0.9045;

CD={'k','b','r','g'}
FalqonCD=0;
%%%%%%%%%%%% Initialization %%%%%%%%%%
HP=HamT(N,1)+HamV(N,1);
HM=HamT(N,1);
HMCD1=HamTVCD(N,1);
%HMCD2=HamNN(N,3,2,1);
%%
CM=(HM*HP-HP*HM);
CMCD1=(HMCD1*HP-HP*HMCD1);
%CMCD2=(HMCD2*HP-HP*HMCD2);
%%
dt=0.01;
T=0:dt:100*dt;
Beta=zeros(length(T),1);
Gamma1=zeros(length(T),1);
%Gamma2=zeros(length(T),1);
Beta(1)=0;
Gamma1(1)=0;
%Gamma2(1)=0;
%% Ground state of T
psiinitialOp=eye(2^N);
for m=-N/4:1:N/4-1
    psiinitial0=zeros(2^N,1);
    for n=1:N;
    psiinitial0=psiinitial0+(exp(-1i*2*pi*(m)*n./N)./sqrt(N)).*Creation(n,N);
    end
    psiinitialOp= psiinitialOp*psiinitial0;
end
vac=[1;0];
for i=1:N-1;
    vac=kron(vac,[1;0]);
end
v=zeros(4^N,length(T));
v(:,1)=kron(psiinitialOp,psiinitialOp)*kron(vac,vac);
energy=zeros(length(T),1);
energy(1)=(v(:,1)'*HP*v(:,1));
%%
%%%%%%%%%%%%%% Evolution %%%%%%%%%%%%%%%

for i=1:length(T)-1
    t=(dt*(i-1));
    h1=Beta(i).*HM;
    h21=Gamma1(i).*HMCD1;
    %h22=Gamma2(i).*HMCD2;
    if FalqonCD==0
       v(:,i+1)=(expm(-1i*dt*h1)*(expm(-1i*dt*HP)*v(:,i)));
       Beta(i+1)=-1i.*(v(:,i)'*CM*v(:,i));
       energy(i+1)=(v(:,i+1)'*HP*v(:,i+1));
    elseif FalqonCD==1
        v(:,i+1)=expm(-1i*dt*h21)*(expm(-1i*dt*h1)*(expm(-1i*dt*HP)*v(:,i)));
        Beta(i+1)=-1i.*(v(:,i)'*CM*v(:,i));
        Gamma1(i+1)=-1i.*(v(:,i)'*CMCD1*v(:,i));
        energy(i+1)=(v(:,i+1)'*HP*v(:,i+1));
    else
        v(:,i+1)=expm(-1i*dt*h22)*(expm(-1i*dt*h1)*(expm(-1i*dt*HP)*v(:,i)));
        Beta(i+1)=-1i.*(v(:,i)'*CM*v(:,i));
        %Gamma2(i+1)=-1i.*(v(:,i)'*CMCD2*v(:,i));
        energy(i+1)=(v(:,i+1)'*HP*v(:,i+1));
    end
end
betadot=zeros(length(Beta),1);
betadot(1)=0;
for i=1:length(Beta)-1
    betadot(i+1)=(Beta(i+1)-Beta(i))./dt;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,3,1)
plot(T,real(energy(:,1)),'Color',CD{FalqonCD+1})
ax1=gca;
xlabel('Time')
ylabel('Energy')
hold on
subplot(1,3,2)
plot(T,real(Beta(:,1)),'Color',CD{FalqonCD+1})
hold on
%plot(T,real(betadot(:,1)),'-o','Color',CD{FalqonCD+1})
ax2=gca;
xlabel('Time')
ylabel('Beta(t)')
hold on
subplot(1,3,3)
if FalqonCD==1
plot(T,real(Gamma1(:,1)),'Color',CD{FalqonCD+1})
%plot(T,real(betadot(:,1)),'Color',CD{FalqonCD+2})
ax3=gca;
xlabel('Time')
ylabel('Gamma(t)')
hold on
else
 %plot(T,real(Gamma2(:,1)),'Color',CD{FalqonCD+1})
%plot(T,real(betadot(:,1)),'Color',CD{FalqonCD+2})
 ax4=gca;
xlabel('Time')
ylabel('Gamma(t)')
 hold on
end








toc