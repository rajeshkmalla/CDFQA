clear all
close all
%% Initialization
N=6;J=-1;hz=-0.4;hx=-0.4;perBC=1;% Hamiltonian parameters
CD={'k',"#0072BD","#D95319","#7E2F8E",'g','c','y'};
pool={'I','Y','ZY','YZ'};
FalqonCD=0;% 0-original FALQON , 1 - FALOQN with Y, 2 - FALOQN with YZ
Colorcount=2;
alpha=1;beta=1;gamma=1;
%%
%% Building Hamiltonian
HP=HamNN(N,3,3,J,perBC)+HamOnsite(N,3,hz)+HamOnsite(N,1,hx);
[EV,D]=eig(HP);
HM=HamOnsite(N,1,1);
HCD0=sparse(eye(2^N));
HCD=HamOnsite(N,2,1);
HCD1=HamNN(N,2,3,1,perBC);
HCD2=HamNN(N,2,1,1,perBC);
%
psi_initial=ones(2^N,1)./sqrt(2^N);
Layers=200;
dt=0.01;

pop0=CDFQA_population_vs_cir_depth(N,HP,HM,HCD0,psi_initial,Layers,dt,alpha);
pop=CDFQA_population_vs_cir_depth(N,HP,HM,HCD,psi_initial,Layers,dt,alpha);
pop1=CDFQA_population_vs_cir_depth(N,HP,HM,HCD1,psi_initial,Layers,dt,alpha);
pop2=CDFQA_population_vs_cir_depth(N,HP,HM,HCD2,psi_initial,Layers,dt,alpha);


sz=size(pop0);
T=length(pop0(1:Layers));
dt=0.01;


D1=diag(D)-D(1,1);
if mod(ceil(max(D1)),2)==0
    ymax=ceil(max(D1));
else
    ymax=ceil(max(D1))+1;
end

%ymax=ceil(max(D1));


P0=zeros((ymax/2)+1,T);
P=zeros((ymax/2)+1,T);
P1=zeros((ymax/2)+1,T);
P2=zeros((ymax/2)+1,T);
n=1;
     for i=1:length(D1)
         D1(i)
    if (2*(n-1)<= D1(i)) && (D1(i)<2*n)
        P0(n,:)=P0(n,:)+pop0(i,:);
        P(n,:)=P(n,:)+pop(i,:);
        P1(n,:)=P1(n,:)+pop1(i,:);
        P2(n,:)=P2(n,:)+pop2(i,:);
        n
    else
        n=n+1
    end

    end

f=figure;
[X,Y] = meshgrid(1:sz(2),0:2/6:ymax/6);
surf(X,Y,P0,'edgecolor','none')

axis([0 200 0 16/6])
view(2)
box on;
colorbar
colormap("turbo")
set(gca, "FontName",'Number of Layers','FontSize',24,'LineWidth',1.5)
xlabel('{\cal L}','Interpreter','latex',FontSize=24)
ylabel('e_P',FontSize=24)
exportgraphics(f,'Fig5/3D-FQA0.jpg')

f=figure;
[X,Y] = meshgrid(1:sz(2),0:2/6:ymax/6);
surf(X,Y,P,'edgecolor','none')

axis([0 200 0 16/6])
view(2)
box on;
colorbar
colormap("turbo")
set(gca, "FontName",'Number of Layers','FontSize',24,'LineWidth',1.5)
xlabel('{\cal L}','Interpreter','latex',FontSize=24)
ylabel('e_P',FontSize=24)
exportgraphics(f,'Fig5/3D-FQA1.jpg')

f=figure;
[X,Y] = meshgrid(1:sz(2),0:2/6:ymax/6);
surf(X,Y,P1,'edgecolor','none')

axis([0 200 0 16/6])
view(2)
box on;
colorbar
colormap("turbo")
set(gca, "FontName",'Number of Layers','FontSize',24,'LineWidth',1.5)
xlabel('{\cal L}','Interpreter','latex',FontSize=24)
ylabel('e_P',FontSize=24)
exportgraphics(f,'Fig5/3D-FQA2.jpg')
f=figure;
[X,Y] = meshgrid(1:sz(2),0:2/6:ymax/6);
surf(X,Y,P2,'edgecolor','none')

axis([0 200 0 16/6])
view(2)
box on;
colorbar
colormap("turbo")
set(gca, "FontName",'Number of Layers','FontSize',24,'LineWidth',1.5)
xlabel('{\cal L}','Interpreter','latex',FontSize=24)
ylabel('e_P',FontSize=24)
exportgraphics(f,'Fig5/3D-FQA3.jpg')