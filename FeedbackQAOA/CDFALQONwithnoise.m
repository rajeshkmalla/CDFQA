clear all
tic
%% Initialization
N=6;J=-1;hz=-0;hx=-0.4;perBC=1;
CD={'k','b','r','g','g','c','y'};
pool={'I','Y','ZY','YZ'};
FalqonCD=3;% 0-original FALQON , 1 - FALOQN with Y, 2 - FALOQN with YZ
Colorcount=2;
alpha=1;beta=1;gamma=1;
%% Building Hamiltonian
HP=HamNN(N,3,3,J,perBC)+HamOnsite(N,3,hz)+HamOnsite(N,1,hx);
%HMCD=HamNN(N,3,2,1,perBC);
HM=HamOnsite(N,1,1);
 %HM=HMCD;

%% Commutators 
CM=(HM*HP-HP*HM);
%% 
dt=0.005;
T=0:dt:2000*dt;
alpha1=T./(40*dt);
Beta=zeros(length(T),1);  
g=zeros(length(T),1); 
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
    h1=Beta(i)*HM;
    HP=HamNN(N,3,3,J,perBC)+HamOnsite(N,3,hz)+HamOnsite(N,1,hx);
    v(:,i+1)=(expm(-1i*dt*beta*h1)*(expm(-1i*dt*alpha*HP)*v(:,i)));
    Beta(i+1)=-1i.*(v(:,i)'*CM*v(:,i));
    energy(i+1)=(v(:,i+1)'*HP*v(:,i+1));
    population(:,i+1)=abs((v(:,i+1)'*EV)').^2;
    end
    %
elseif FalqonCD==1
    HMCD=HamOnsite(N,2,1);
     CMCD=(HMCD*(HP))-((HP)*HMCD);
    %
    for i=1:length(T)-1
    t=(dt*(i-1));
    h1=Beta(i).*HM;
    h2=Gamma(i).*HMCD;
    %h2=(Gamma(i)+abs(energy(i))./2).*HMCD;
    % hz=(-0.5+rand(1))*(1/i.^2);
    HP=HamNN(N,3,3,J,perBC)+HamOnsite(N,3,hz)+HamOnsite(N,1,hx);
    v(:,i+1)=expm(-1i*dt*gamma*h2)*(expm(-1i*dt*beta*h1)*(expm(-1i*dt*(HP+HamOnsite(N,3,hz)))*v(:,i)));
    Beta(i+1)=-1i.*(v(:,i)'*CM*v(:,i));
    Gamma(i+1)=-1i.*(v(:,i)'*CMCD*v(:,i));
    energy(i+1)=(v(:,i+1)'*HP*v(:,i+1));
    population(:,i+1)=abs((v(:,i+1)'*EV)').^2;
    end
    %
elseif FalqonCD==2
    HMCD=HamNN(N,2,3,1,perBC);
    CMCD=(HMCD*(HP))-((HP)*HMCD);
    %
    for i=1:(length(T))-1
    t=(dt*(i-1));
    %CMCD=(HMCD*(HP+Beta(i).*HM))-((HP+Beta(i).*HM)*HMCD);
    % if abs(Beta(i))<1
    %     Beta(i)=sign(real(Beta(i)))*sqrt(abs(Beta(i)));
    % end
    % if abs(Gamma(i))<1
    %     Gamma(i)=sign(real(Gamma(i)))*sqrt(abs(Gamma(i)));
    % end
     h1=Beta(i).*HM;
     %h1=(Beta(i)+abs(energy(i)-D(1,1))./2).*HM;
     %g(i)= (energy(i)-D(1,1))/2;
     h2=(Gamma(i)).*HMCD;
    %h2=(Gamma(i)+abs(energy(i)-D(1,1))./2).*HMCD;
    %hz=(-0.5+rand(1))*(1/i.^2);
    %HP=HamNN(N,3,3,J,perBC)+HamOnsite(N,3,hz)+HamOnsite(N,1,hx);
    % if(i>=40)
    %     alpha=0;
    % end
    v(:,i+1)=expm(-1i*dt*gamma*h2)*(expm(-1i*dt*beta*h1)*(expm(-1i*dt*(alpha*HP))*v(:,i)));
    Beta(i+1)=-1i.*(v(:,i)'*CM*v(:,i));
    Gamma(i+1)=-1i.*(v(:,i)'*CMCD*v(:,i));
    energy(i+1)=(v(:,i+1)'*HP*v(:,i+1));
    population(:,i+1)=abs((v(:,i+1)'*EV)').^2;
    end
    %
 elseif FalqonCD==3
    HMCD=HamNN(N,2,2,1,perBC)+HamOnsite(N,2,1);
    CMCD=(HMCD*(HP)-(HP)*HMCD);
    %
    for i=1:length(T)-1
    t=(dt*(i-1));
    h1=Beta(i).*HM;
    h2=Gamma(i).*HMCD;
    v(:,i+1)=expm(-1i*dt*gamma*h2)*(expm(-1i*dt*beta*h1)*(expm(-1i*dt*alpha*HP)*v(:,i)));
    Beta(i+1)=-1i.*(v(:,i)'*CM*v(:,i))*beta;
    Gamma(i+1)=-1i.*(v(:,i)'*CMCD*v(:,i));
    energy(i+1)=(v(:,i+1)'*HP*v(:,i+1));
    population(:,i+1)=abs((v(:,i+1)'*EV)').^2;
    end
else
    disp('Invalid CD option')
end 






%%%%%%%%%%%%%%%%%%NEW run with constant 
% energyconst=abs(energy(i));
% for i=length(T)+1:2*length(T)-1
%     t=(dt*(i-1));
%     h1=sign(Beta(i))*energyconst.*HM;
%     v(:,i+1)=(expm(-1i*dt*beta*h1)*(expm(-1i*dt*alpha*HP)*v(:,i)));
%     Beta(i+1)=-1i.*(v(:,i)'*CM*v(:,i));
%     energy(i+1)=(v(:,i+1)'*HP*v(:,i+1));
%     population(:,i+1)=abs((v(:,i+1)'*EV)').^2;
% end


% Plots
subplot(1,3,1)
plot(T,real(energy(:,1))-D(1,1),'Color',CD{FalqonCD+Colorcount+1},'LineWidth',2)
hold on
% for i=1:2^6
% plot(T,real(energy(:,i))-D(1,1),'Color',CD{FalqonCD+Colorcount+1},'LineWidth',2)
% hold on
% end
ax1=gca;
xlabel('Time')
ylabel('Energy')
hold on
subplot(1,3,2)
plot(T,real(Beta(:,1)),'Color',CD{FalqonCD+Colorcount+1},'LineWidth',2)
hold on
%plot(T,real(betadot(:,1)),'-o','Color',CD{FalqonCD+1})
ax2=gca;
xlabel('Time')
ylabel('Beta(t)')
hold on
subplot(1,3,3)
plot(T,real(Gamma(:,1)),'Color',CD{FalqonCD+Colorcount+1},'LineWidth',2)
hold on
%plot(T,real(betadot(:,1)),'Color',CD{FalqonCD+2})
ax3=gca;
xlabel('Time')
ylabel('Gamma(t)')

toc
ylabel('Energy')
hold on
subplot(1,3,2)
plot(T,real(Beta(:,1)),'Color',CD{FalqonCD+Colorcount+1},'LineWidth',2)
hold on
%plot(T,real(betadot(:,1)),'-o','Color',CD{FalqonCD+1})
ax2=gca;
xlabel('Time')
ylabel('Beta(t)')
hold on
subplot(1,3,3)
plot(T,real(Gamma(:,1)),'Color',CD{FalqonCD+Colorcount+1},'LineWidth',2)
hold on
%plot(T,real(betadot(:,1)),'Color',CD{FalqonCD+2})
ax3=gca;
xlabel('Time')
ylabel('Gamma(t)')


% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Initialization
% N=6;J=-1;hz=-0;hx=+0;perBC=1;
% CD={'k','b','r','g','m','c','y'};
% pool={'I','Y','ZY','YZ'};
% FalqonCD=1%oginal FALQON , 1 - FALOQN with Y, 2 - FALOQN with YZ
% Colorcount=1;
% alpha=1;beta=1;gamma=1;
% %% Building Hamiltonian
% HP=HamNN(N,3,3,J,perBC)+HamOnsite(N,3,hz)+HamOnsite(N,1,hx);
% %HMCD=HamNN(N,3,2,1,perBC);
% HM=HamOnsite(N,1,1);
%  %HM=HMCD;
% 
% %% Commutators 
% CM=(HM*HP-HP*HM);
% %% 
% dt=0.01;
% T=0:dt:40*dt;
% Beta=zeros(length(T),1);
% Gamma=zeros(length(T),1);
% Beta(1)=0;
% Gamma(1)=0;
% 
% 
% %% Population initialization
% v1=zeros(2^N,length(T));
% population=zeros(2^N,length(T));
% energy1=zeros(length(T),1);
% v1(:,1)=v(:,length(T));
% energy1(1)=(v1(:,1)'*HP*v1(:,1));
% [EV,D]=eig(HP);
% population(:,1)=abs((v1(:,1)'*EV)').^2;

% %% Time evolution
% %%%%%%%%%%%%%%%%%%NEW run with constant 
% energyconst=abs(energy1(1));
% Beta(1)=Beta(length(T))
% for i=1:length(T)-1
%     t=(dt*(i-1));
%     h1=sign(Beta(i))*5*HM;
%     v1(:,i+1)=(expm(-1i*dt*beta*h1)*(expm(-1i*dt*alpha*HP)*v1(:,i)));
%     Beta(i+1)=-1i.*(v1(:,i)'*CM*v1(:,i));
%     energy1(i+1)=(v1(:,i+1)'*HP*v1(:,i+1));
%     population(:,i+1)=abs((v1(:,i+1)'*EV)').^2;
% end
% 
% subplot(1,3,1)
% plot(1:length(T),real(energy1(:,1))-D(1,1),'Color',CD{FalqonCD+Colorcount+1},'LineWidth',2)
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
% 
% toc
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
% 
toc