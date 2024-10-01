clear all

%% Initialization
N=6;J=-1;hz=-0;hx=-0;perBC=1;
CD={'k','b','r','g','m','c','y'};
FalqonCD=0;% 0 - Original FALQON , 1 - FALOQN with Y, 2 - FALOQN with YZ
Colorcount=1;
%% Building Hamiltonian
HP=HamNN(N,3,3,J,perBC)+HamOnsite(N,3,hz)+HamOnsite(N,1,hx);
HM=HamOnsite(N,1,1);


%% Commutators 
CM=(HM*HP-HP*HM);
%% 
dt=0.001;
T=0:dt:200*dt;
Beta=zeros(length(T),1);
Gamma=zeros(length(T),1);
Beta(1)=0;
Gamma(1)=0;

%% Population initialization
v=zeros(2^N,length(T));
population=zeros(2^N,length(T));
v(:,1)=ones(2^N,1)./sqrt(2^N);
energy=zeros(length(T),1);
energyderivative=zeros(length(T)-1,1);
energy(1)=(v(:,1)'*HP*v(:,1));
[EV,D]=eig(HP);
population(:,1)=abs((v(:,1)'*EV)').^2;

%% Time evolution
%gammajump=5;
FalqonCD=2;
HMCD=HamNN(N,3,2,1,perBC);
CMCD=(HMCD*HP-HP*HMCD);
count=1;
secondcount=1;
for i=1:length(T)-1
    t=(dt*(i-1));
    h1=Beta(i).*HM;
    h2=Gamma(i).*HMCD;
    if FalqonCD==2
    v(:,i+1)=expm(-1i*dt*h2)*(expm(-1i*dt*h1)*(expm(-1i*dt*HP)*v(:,i)));
    Beta(i+1)=-1i.*(v(:,i)'*CM*v(:,i));    
    if count==2
        count;
        gammafix;
        if secondcount==1
        sign(real(Gamma(i)))
        secondcount=secondcount+1;
        end
        Gamma(i+1)=sign(real(Gamma(i)))*abs(gammafix);
        %Beta(i+1)=sign(Beta(i))*; 
         %Gamma(i+1)=-1i.*(v(:,i)'*CMCD*v(:,i));
    else
    Gamma(i+1)=-1i.*(v(:,i)'*CMCD*v(:,i));
    end

    elseif FalqonCD==0
    v(:,i+1)=(expm(-1i*dt*h1)*(expm(-1i*dt*HP)*v(:,i)));
    Beta(i+1)=-1i.*(v(:,i)'*CM*v(:,i));
    Beta(i+1)=sign(Beta(i+1))*1;
    Gamma(i+1)=0;
    end
    
    
    % if count==1
    %         Beta(i+1)=sign(Beta(i+1))*10;
    % end
    % count=count+1;
    % end
     energy(i+1)=(v(:,i+1)'*HP*v(:,i+1));
     energyderivative(i)=abs(energy(i+1)-energy(i));
    
     if abs(Gamma(i+1))< 0.1 && i>2 && count==1
         gammafix=abs(Gamma(i));
         sign(real(Gamma(i+1)))
         sign(real(Gamma(i+1)))*abs(gammafix)
         Gamma(i+1)=sign(real(Gamma(i+1)))*abs(gammafix);
         count=count+1;
         FalqonCD=0;
        % v(:,i+1)=(expm(-1i*5*dt*HP)*v(:,i+1));
       %Beta(i+1)=sign(Beta(i+1))*1;
       %Gamma(i+1)=sign(Gamma(i+1))*gammafix;
     end

     
     % if abs(energy(i+1)-energy(i))<=0.1
    %    FalqonCD=0;
    %   % Gamma(i+1)=gammajump*sign(Gamma(i+1));
    % else
    %     FalqonCD=2;
    % end
end



% if FalqonCD==0
%     disp('Original Falqon')
%     %
%     for i=1:length(T)-1
%     t=(dt*(i-1));
%     h1=Beta(i).*HM;
%     v(:,i+1)=(expm(-1i*dt*h1)*(expm(-1i*dt*HP)*v(:,i)));
%     Beta(i+1)=-1i.*(v(:,i)'*CM*v(:,i));
%     energy(i+1)=(v(:,i+1)'*HP*v(:,i+1));
%     population(:,i+1)=abs((v(:,i+1)'*EV)').^2;
%     end
%     %
% elseif FalqonCD==1
%     HMCD=HamOnsite(N,2,1);
%     CMCD=(HMCD*HP-HP*HMCD);
%     %
%     for i=1:length(T)-1
%     t=(dt*(i-1));
%     h1=Beta(i).*HM;
%     h2=Gamma(i).*HMCD;
%     v(:,i+1)=expm(-1i*dt*h2)*(expm(-1i*dt*h1)*(expm(-1i*dt*HP)*v(:,i)));
%     Beta(i+1)=-1i.*(v(:,i)'*CM*v(:,i));
%     Gamma(i+1)=-1i.*(v(:,i)'*CMCD*v(:,i));
%     energy(i+1)=(v(:,i+1)'*HP*v(:,i+1));
%     population(:,i+1)=abs((v(:,i+1)'*EV)').^2;
%     end
%     %
% elseif FalqonCD==2
%     HMCD=HamNN(N,3,2,1,perBC);
%     CMCD=(HMCD*HP-HP*HMCD);
%     %
%     for i=1:length(T)-1
%     t=(dt*(i-1));
%     h1=Beta(i).*HM;
%     h2=Gamma(i).*HMCD;
%     v(:,i+1)=expm(-1i*dt*h2)*(expm(-1i*dt*h1)*(expm(-1i*dt*HP)*v(:,i)));
%     Beta(i+1)=-1i.*(v(:,i)'*CM*v(:,i));
%     Gamma(i+1)=-1i.*(v(:,i)'*CMCD*v(:,i));
%     energy(i+1)=(v(:,i+1)'*HP*v(:,i+1));
%     population(:,i+1)=abs((v(:,i+1)'*EV)').^2;
%     end
%     %
%  elseif FalqonCD==3
%     HMCD=HamNN(N,2,3,1,perBC);
%     CMCD=(HMCD*HP-HP*HMCD);
%     %
%     for i=1:length(T)-1
%     t=(dt*(i-1));
%     h1=Beta(i).*HM;
%     h2=Gamma(i).*HMCD;
%     v(:,i+1)=expm(-1i*dt*h2)*(expm(-1i*dt*h1)*(expm(-1i*dt*HP)*v(:,i)));
%     Beta(i+1)=-1i.*(v(:,i)'*CM*v(:,i));
%     Gamma(i+1)=-1i.*(v(:,i)'*CMCD*v(:,i));
%     energy(i+1)=(v(:,i+1)'*HP*v(:,i+1));
%     population(:,i+1)=abs((v(:,i+1)'*EV)').^2;
%     end
% else
%     disp('Invalid CD option')
% end 

% Plots
subplot(1,3,1)
plot(T,real(energy(:,1))-D(1,1),'Color',CD{FalqonCD+Colorcount+1},'LineWidth',2)
hold on
%plot(T(2:length(T)),energyderivative(:),'Color',CD{FalqonCD+Colorcount+1},'LineWidth',2)

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

% %% Population plots
% videofile = VideoWriter('FALQON-0','MPEG-4');
% open(videofile);
% fig1=figure(1);
% for i=1:length(T)
% plot(1:2^N,population(:,i))
% hold on
% axis([0 2^N 0 1])
% time = (i-1)*dt;
% text(2^(N/2),0.7,num2str(time))
% %pause(0.1)
% F=getframe(fig1);
% hold off
% writeVideo(videofile,F);
% end
% close(videofile);
