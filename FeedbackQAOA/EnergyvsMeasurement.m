clear all
close all
%% Initialization
N=6;J=-1;hz=-0.4;hx=-0;perBC=1;% Hamiltonian parameters
CD={'k',"#0072BD","#D95319","#7E2F8E",'g','c','y'};
pool={'I','Y','ZY','YZ'};
FalqonCD=3;% 0-original FALQON , 1 - FALOQN with Y, 2 - FALOQN with YZ
Colorcount=2;
alpha=1;beta=1;gamma=1;
%%
%% Building Hamiltonian
HP=HamNN(N,3,3,J,perBC)+HamOnsite(N,3,hz)+HamOnsite(N,1,hx);
%Hpert=HamOnsite(N,3,hz);
%HMCD=HamNN(N,3,2,1,perBC);
HM=HamOnsite(N,1,1);
 %HM=HMCD;

%% Commutators 
CM=(HM*HP-HP*HM);
%% 
dt=0.01;
T=0:dt:800*dt;
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
f=figure;
for FalqonCD=0:3
    if FalqonCD==0
    disp('Original Falqon')
    %
    for i=1:length(T)-1
    t=(dt*(i-1));
    h1=Beta(i)*HM;
    HP=HamNN(N,3,3,J,perBC)+HamOnsite(N,3,hz)+HamOnsite(N,1,hx);
    v(:,i+1)=(fastExpm(-1i*dt*beta*h1)*(fastExpm(-1i*dt*alpha*HP)*v(:,i)));
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
    HP=HamNN(N,3,3,J,perBC)+HamOnsite(N,3,hz)+HamOnsite(N,1,hx);
    v(:,i+1)=fastExpm(-1i*dt*gamma*h2)*(fastExpm(-1i*dt*beta*h1)*(fastExpm(-1i*dt*(HP))*v(:,i)));
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
     h1=Beta(i).*HM;
     h2=(Gamma(i)).*HMCD;
    v(:,i+1)=fastExpm(-1i*dt*gamma*h2)*(fastExpm(-1i*dt*beta*h1)*(fastExpm(-1i*dt*(alpha*HP))*v(:,i)));
    Beta(i+1)=-1i.*(v(:,i)'*CM*v(:,i));
    Gamma(i+1)=-1i.*(v(:,i)'*CMCD*v(:,i));
    energy(i+1)=(v(:,i+1)'*HP*v(:,i+1));
    population(:,i+1)=abs((v(:,i+1)'*EV)').^2;
    end
    %
 elseif FalqonCD==3
    HMCD=HamNN(N,2,1,1,perBC);
    CMCD=(HMCD*(HP)-(HP)*HMCD);
    %
    for i=1:length(T)-1
    t=(dt*(i-1));
    h1=Beta(i).*HM;
    h2=Gamma(i).*HMCD;
    v(:,i+1)=fastExpm(-1i*dt*gamma*h2)*(fastExpm(-1i*dt*beta*h1)*(fastExpm(-1i*dt*alpha*HP)*v(:,i)));
    Beta(i+1)=-1i.*(v(:,i)'*CM*v(:,i))*beta;
    Gamma(i+1)=-1i.*(v(:,i)'*CMCD*v(:,i));
    energy(i+1)=(v(:,i+1)'*HP*v(:,i+1));
    population(:,i+1)=abs((v(:,i+1)'*EV)').^2;
    end
else
    disp('Invalid CD option')
    end
%%
mdata(:,FalqonCD+1)=(real(energy(:,1))-D(1,1))/N;
% plot(1:length(T),(real(energy(:,1))-D(1,1))/N,'-o','Color',CD{FalqonCD+1},'LineWidth',2)
% hold on
%  %plot(1:length(T),Gamma(:),'o-','Color',CD{FalqonCD+1},'LineWidth',2)
%  %hold on
% % plot(1:length(T),Gamma(:),'x','Color',CD{FalqonCD+1},'LineWidth',2)
% % hold on
% %scatter(T,real(energy(:,1)),25,Marker="o",MarkerEdgeColor=CD{FalqonCD+1},LineWidth=2)
% grid on
% legend('I','Y','YZ','YX')
% %legend('I,Y','YZ')
% set(gca, "FontName",'Number of Layers','FontSize',24,'LineWidth',1.5)
% xlabel('{\cal L}','Interpreter','latex',FontSize=28)
% ylabel('e_P',FontSize=28)
% %ylabel('\beta(t)',FontWeight='bold',FontSize=24)
% %ylabel('\gamma(t)',FontWeight='bold',FontSize=24)


end
Xdata=zeros(length(T)*8,4);
for t=1:length(T)
for i=1:4
if i==1
Xdata(2*t,i)=mdata(t,i);
end

if i==2
Xdata(4*t,i)=mdata(t,i);
end

if i==3
Xdata(4*t,i)=mdata(t,i);
end

if i==4
Xdata(8*t,i)=mdata(t,i);
end
end
end

%Xdata1=nonzeros(Xdata);

    plot(1:length(T)/2,Xdata(1:length(T)/2,1),'o','Color',CD{1},'LineWidth',2)
    hold on
    plot(1:length(T),Xdata(1:length(T),2),'o','Color',CD{2},'LineWidth',2)
    plot(1:length(T),Xdata(1:length(T),3),'o','Color',CD{3},'LineWidth',2)
    plot(1:length(T)*2,Xdata(1:length(T)*2,4),'o','Color',CD{4},'LineWidth',2)
    
    h1=plot(length(T)/2:length(T)*2,Xdata(length(T)/2:length(T)*2,1),'o','Color',CD{1},'LineWidth',2)
    hold on
    h2=plot(length(T):length(T)*2,Xdata(length(T):length(T)*2,2),'o','Color',CD{2},'LineWidth',2)
    h3=plot(length(T):length(T)*2,Xdata(length(T):length(T)*2,3),'o','Color',CD{3},'LineWidth',2)
set(h1,'Color', 1-0.5*(1-get(h1,'Color')))
set(h2,'Color', 1-0.5*(1-get(h2,'Color')))
set(h3,'Color', 1-0.5*(1-get(h3,'Color')))
grid on
legend('I(\times 2)','Y(\times 4)','YZ(\times 4)','YX(\times 8)')
%legend('I,Y','YZ')
set(gca, "FontName",'Number of Layers','FontSize',24,'LineWidth',1.5)
%xlabel('n_m','Interpreter','latex',FontSize=28)
xlabel('N_{meas}',FontSize=28)
ylabel('e_P',FontSize=28)
axis([0.01 length(T)*2 0.01 1.5])
%exportgraphics(f,'TFI114.jpg')
