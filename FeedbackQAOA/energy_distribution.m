clear all
tic
%% Initialization
N=6;J=-1;hz=-0.4;hx=-0.4;perBC=1;
CD={'k','b','r','g','m','c','y'};
pool={'I','Y','ZY','YZ'};
FalqonCD=3% 0 - Original FALQON , 1 - FALOQN with Y, 2 - FALOQN with YZ
Colorcount=1;
beta=1;gamma=1;alpha=1;
%% Building Hamiltonian
HP=HamNN(N,3,3,J,perBC)+HamOnsite(N,3,hz)+HamOnsite(N,1,hx);
%HMCD=HamNN(N,3,2,1,perBC);
HM=HamOnsite(N,1,1);
 %HM=HMCD;

%% Commutators 
CM=(HM*HP-HP*HM);
%% 
dt=0.01;
T=0:dt:200*dt;
Beta=zeros(length(T),1);
Gamma=zeros(length(T),1);
Beta(1)=0;
Gamma(1)=0;


%% Population initialization
v=zeros(2^N,length(T));
population=zeros(2^N,length(T));
energy=zeros(length(T),1);
v(:,1)=ones(2^N,1)./sqrt(2^N);
energy(1)=(v(:,1)'*HP*v(:,1));
[EV,D]=eig(HP);
population(:,1)=abs((v(:,1)'*EV)').^2;

%% Time evolution
if FalqonCD==0
    disp('Original Falqon')
    %
    for i=1:length(T)-1
    t=(dt*(i-1));
    h1=Beta(i).*HM;
    v(:,i+1)=(expm(-1i*dt*h1)*(expm(-1i*dt*alpha*HP)*v(:,i)));
    Beta(i+1)=-1i.*(v(:,i)'*CM*v(:,i));
    energy(i+1)=(v(:,i+1)'*HP*v(:,i+1));
    population(:,i+1)=abs((v(:,i+1)'*EV)').^2;
    end
    %
elseif FalqonCD==1
    HMCD=HamOnsite(N,2,1);
    CMCD=(HMCD*HP-HP*HMCD);
    %
    for i=1:length(T)-1
    t=(dt*(i-1));
    h1=Beta(i).*HM;
    h2=Gamma(i).*HMCD;
    v(:,i+1)=expm(-1i*dt*h2)*(expm(-1i*dt*h1)*(expm(-1i*dt*alpha*HP)*v(:,i)));
    Beta(i+1)=-1i.*(v(:,i)'*CM*v(:,i));
    Gamma(i+1)=-1i.*(v(:,i)'*CMCD*v(:,i));
    energy(i+1)=(v(:,i+1)'*HP*v(:,i+1));
    population(:,i+1)=abs((v(:,i+1)'*EV)').^2;
    end
    %
elseif FalqonCD==2
    HMCD=HamNN(N,2,3,1,perBC);
    CMCD=(HMCD*HP-HP*HMCD);
    %
    for i=1:length(T)-1
    t=(dt*(i-1));
    % if abs(Beta(i))<1
    %     Beta(i)=sign(real(Beta(i)))*sqrt(abs(Beta(i)));
    % end
    % if abs(Gamma(i))<1
    %     Gamma(i)=sign(real(Gamma(i)))*sqrt(abs(Gamma(i)));
    % end
     h1=Beta(i).*HM;
    h2=Gamma(i).*HMCD;
    v(:,i+1)=expm(-1i*dt*h2)*(expm(-1i*dt*h1)*(expm(-1i*dt*alpha*HP)*v(:,i)));
    Beta(i+1)=-1i.*(v(:,i)'*CM*v(:,i));
    Gamma(i+1)=-1i.*(v(:,i)'*CMCD*v(:,i));
    energy(i+1)=(v(:,i+1)'*HP*v(:,i+1));
    population(:,i+1)=abs((v(:,i+1)'*EV)').^2;
    end
    %
 elseif FalqonCD==3
    HMCD=HamNN(N,2,1,1,perBC);
    CMCD=(HMCD*HP-HP*HMCD);
    %
    for i=1:length(T)-1
    t=(dt*(i-1));
    h1=Beta(i).*HM;
    h2=Gamma(i).*HMCD;
    v(:,i+1)=expm(-1i*dt*h2)*(expm(-1i*dt*h1)*(expm(-1i*dt*HP)*v(:,i)));
    Beta(i+1)=-1i.*(v(:,i)'*CM*v(:,i))*beta;
    Gamma(i+1)=-1i.*(v(:,i)'*CMCD*v(:,i));
    energy(i+1)=(v(:,i+1)'*HP*v(:,i+1));
    population(:,i+1)=abs((v(:,i+1)'*EV)').^2;
    end
else
    disp('Invalid CD option')
end 


% fenergy= sprintf('dataCD/MFIM/FalqonCDenergy%d.%f.%f.%f.%d.mat',N,J,hz,hx,FalqonCD);
% save(fenergy,'energy')
% 
% fBeta= sprintf('dataCD/MFIM/FalqonCDBeta%d.%f.%f.%f.%d.mat',N,J,hz,hx,FalqonCD);
% save(fBeta,'Beta')
% 
% fGamma= sprintf('dataCD/MFIM/FalqonCDGamma%d.%f.%f.%f.%d.mat',N,J,hz,hx,FalqonCD);
% save(fGamma,'Gamma')
% 
 fpopulation= sprintf('dataCD/FalqonCDpopulation%d.%f.%f.%f.%d.mat',N,J,hz,hx,FalqonCD);
 save(fpopulation,'population')



% Plots
%subplot(1,3,1)
% scatte(1:length(T),real(energy(:,1))-D(1,1),'Color',CD{FalqonCD+Colorcount+1},'LineWidth',2)
% hold on
% ax1=gca;
% xlabel('Time')
% ylabel('Energy')
% hold on
% subplot(1,3,2)
% plot(T,real(Beta(:,1)),'Color',CD{FalqonCD+Colorcount+1},'LineWidth',2)
% hold on
% %plot(T,real(betadot(:,1)),'-o','Color',CD{FalqonCD+1})
% ax2=gca;
% xlabel('Time')
% ylabel('Beta(t)')
% hold on
% subplot(1,3,3)
% plot(T,real(Gamma(:,1)),'Color',CD{FalqonCD+Colorcount+1},'LineWidth',2)
% hold on
% %plot(T,real(betadot(:,1)),'Color',CD{FalqonCD+2})
% ax3=gca;
% xlabel('Time')
% ylabel('Gamma(t)')

toc