clear all
close all
tic
CD1={'','o-'}
%% Initialization
f=figure;
count=0;
% for hx=-0.4:0.4:0
%     %CD1{count}=sprintf('N=%d',N);count=count+1;
% end
for hx=-0.4:0.4:0
    count=count+1;
    N=4;
J=-1;hz=-0.4;perBC=1;
CD={'k',"#0072BD","#D95319","#7E2F8E",'g','c','c','y'};
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
T=0:dt:100*dt;
alpha1=T./(40*dt);
Beta=zeros(length(T),1);  
g=zeros(length(T),1); 
Gamma=zeros(length(T),1);
Beta(1)=0;
Gamma(1)=0;


%% Population initialization
energy=zeros(length(T),1);
v(:,1)=sparse(ones(2^N,1)./sqrt(2^N));
energy(1)=(v(:,1)'*HP*v(:,1));
D=eigs(HP);
E0=min(D);
%population(:,1)=abs((v(:,1)'*EV)').^2;


for FalqonCD=0:3
    if FalqonCD==0
    disp('Original Falqon')
    %
    for i=1:length(T)-1
        
    t=(dt*(i-1));
    h1=Beta(i)*HM;
    %HP=HamNN(N,3,3,J,perBC)+HamOnsite(N,3,hz)+HamOnsite(N,1,hx);
    v(:,i+1)=fastExpm(-1i*dt*alpha*HP)*v(:,i);
    v(:,i+1)=fastExpm(-1i*dt*beta*h1)*v(:,i+1);
    Beta(i+1)=-1i.*(v(:,i)'*CM*v(:,i));
    energy(i+1)=(v(:,i+1)'*HP*v(:,i+1));
    %population(:,i+1)=abs((v(:,i+1)'*EV)').^2;
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
    %dd*9
    % 'x48P=HamNN(N,3,3,J,perBC)+HamOnsite(N,3,hz)+HamOnsite(N,1,hx);
    v(:,i+1)=fastExpm(-1i*dt*(HP))*v(:,i);
    v(:,i+1)=fastExpm(-1i*dt*beta*h1)*v(:,i+1);
    v(:,i+1)=fastExpm(-1i*dt*gamma*h2)*v(:,i+1);
    Beta(i+1)=-1i.*(v(:,i)'*CM*v(:,i));
    Gamma(i+1)=-1i.*(v(:,i)'*CMCD*v(:,i));
    energy(i+1)=(v(:,i+1)'*HP*v(:,i+1));
    %population(:,i+1)=abs((v(:,i+1)'*EV)').^2;
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
    v(:,i+1)=fastExpm(-1i*dt*(HP))*v(:,i);
    v(:,i+1)=fastExpm(-1i*dt*beta*h1)*v(:,i+1);
    v(:,i+1)=fastExpm(-1i*dt*gamma*h2)*v(:,i+1);
    Beta(i+1)=-1i.*(v(:,i)'*CM*v(:,i));
    Gamma(i+1)=-1i.*(v(:,i)'*CMCD*v(:,i));
    energy(i+1)=(v(:,i+1)'*HP*v(:,i+1));
    %population(:,i+1)=abs((v(:,i+1)'*EV)').^2;
    end
    %
 elseif FalqonCD==3
    HMCD=HamNN(N,2,1,1,perBC);%+ HamNN(N,2,1,1,perBC);
    CMCD=(HMCD*(HP)-(HP)*HMCD);
    %
    for i=1:length(T)-1
    t=(dt*(i-1));
    h1=Beta(i).*HM;
    h2=Gamma(i).*HMCD;
    v(:,i+1)=fastExpm(-1i*dt*(HP))*v(:,i);
    v(:,i+1)=fastExpm(-1i*dt*beta*h1)*v(:,i+1);
    v(:,i+1)=fastExpm(-1i*dt*gamma*h2)*v(:,i+1);    
    Beta(i+1)=-1i.*(v(:,i)'*CM*v(:,i))*beta;
    Gamma(i+1)=-1i.*(v(:,i)'*CMCD*v(:,i));
    energy(i+1)=(v(:,i+1)'*HP*v(:,i+1));
    %population(1,i+1)=abs((v(1,i+1)'*EV)').^2;
    end
else
    disp('Invalid CD option')
    end
%%

%subplot(3,1,1)
plot(1:length(T),real(energy(:,1)-E0)/N,CD1{count},'Color',CD{FalqonCD+1},'LineWidth',3)
hold on
% subplot(3,1,2)
%  plot(1:length(T),Beta(:),'o-','Color',CD{FalqonCD+1},'LineWidth',2)
% hold on
% subplot(3,1,3)
% plot(1:length(T),Gamma(:),'x','Color',CD{FalqonCD+1},'LineWidth',2)
% hold on
%scatter(T,real(energy(:,1)),25,Marker="o",MarkerEdgeColor=CD{FalqonCD+1},LineWidth=2)
grid on
%legend(CD1{:})
set(gca, "FontName",'Number of Layers','FontSize',24,'LineWidth',1.5)
xlabel('Number of Layers',FontSize=24)
ylabel('\langle\Delta E\rangle_t/N',FontSize=24)
%ylabel('\beta(t)',FontWeight='bold',FontSize=24)
%ylabel('\gamma(t)',FontWeight='bold',FontSize=24)

end
clear v
end
toc
%exportgraphics(f,'Plots/LFIMFI.jpg')
