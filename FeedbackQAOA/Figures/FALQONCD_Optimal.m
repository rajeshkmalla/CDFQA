clear all
%close all
tic
N=8;
J=-1;
hz=-0.5;
hx=-0.0;
CD={'k','b','r','g'}
FalqonCD=2;
%%%%%%%%%%%% Initialization %%%%%%%%%%
HP=HamNN(N,3,3,J)+HamOnsite(N,3,hz)+HamOnsite(N,1,hx);
HM=HamOnsite(N,1,1);
HMCD1=HamOnsite(N,2,1);
HMCD2=HamNN(N,3,2,1);
%%
CM=(HM*HP-HP*HM);
CMCD1=(HMCD1*HP-HP*HMCD1);
CMCD2=(HMCD2*HP-HP*HMCD2);
%%
dt=0.01;
T=0:dt:100*dt;
Beta=zeros(length(T),1);
Gamma1=zeros(length(T),1);
Gamma2=zeros(length(T),1);
Beta(1)=0;
Gamma1(1)=0;
Gamma2(1)=0;
%%
v=zeros(2^N,length(T));
v(:,1)=ones(2^N,1)./sqrt(2^N);
energy=zeros(length(T),1);
energy(1)=(v(:,1)'*HP*v(:,1));
%%
%%%%%%%%%%%%%% Evolution %%%%%%%%%%%%%%%

for i=1:length(T)-1
    t=(dt*(i-1));
    h1=Beta(i).*HM;
    h21=Gamma1(i).*HMCD1;
    h22=Gamma2(i).*HMCD2;
    if FalqonCD==0
       v(:,i+1)=(expm(-1i*dt*h1)*(expm(-1i*dt*HP)*v(:,i)));
       Beta(i+1)=-1i.*(v(:,i)'*CM*v(:,i));
       energy(i+1)=(v(:,i+1)'*HP*v(:,i+1));
    elseif FalqonCD==1
        v(:,i+1)=expm(-1i*dt*h21)*(expm(-1i*dt*h1)*(expm(-1i*dt*HP)*v(:,i)));
        T1=real(1i.*(v(:,i)'*CM*v(:,i)));
        T2=real(1i.*(v(:,i)'*CMCD1*v(:,i)));
        Beta(i+1)=-;
        Gamma1(i+1)=-1i.*(v(:,i)'*CMCD1*v(:,i));
        energy(i+1)=(v(:,i+1)'*HP*v(:,i+1));
    else
        v(:,i+1)=expm(-1i*dt*h22)*(expm(-1i*dt*h1)*(expm(-1i*dt*HP)*v(:,i)));
        Beta(i+1)=-1i.*(v(:,i)'*CM*v(:,i));
        Gamma2(i+1)=-1i.*(v(:,i)'*CMCD2*v(:,i));
        energy(i+1)=(v(:,i+1)'*HP*v(:,i+1));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,1,1)
plot(T,real(energy(:,1)),'Color',CD{FalqonCD+1})
hold on
subplot(3,1,2)
plot(T,real(Beta(:,1)),'Color',CD{FalqonCD+1})
hold on
subplot(3,1,3)
if FalqonCD==1
plot(T,real(Gamma1(:,1)),'Color',CD{FalqonCD+1})
hold on
else
 plot(T,real(Gamma2(:,1)),'Color',CD{FalqonCD+1})
 hold on
end








toc