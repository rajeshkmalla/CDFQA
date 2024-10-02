clear all
close all
%% Initialization
N=6;J=-1;hz=-0;hx=-0.4;perBC=1;% Hamiltonian parameters
CDcolor={'k',"#0072BD","#D95319","#7E2F8E",'g','c','y'};
pool={'I','Y','ZY','YZ'};
FCD=0;% 0-original FALQON , 1 - FALOQN with Y, 2 - FALOQN with YZ

alpha=1;
%%
%% Building Hamiltonian
HP=HamNN(N,3,3,J,perBC)+HamOnsite(N,3,hz)+HamOnsite(N,1,hx);
[EV,D]=eig(HP);
HM=HamOnsite(N,1,1);

HCD0=sparse(eye(2^N));
HCD1=HamOnsite(N,2,1);
HCD2=HamNN(N,2,3,1,perBC);
HCD3=HamNN(N,2,1,1,perBC);



%
psi_initial=ones(2^N,1)./sqrt(2^N);
Layers=200;
dt=0.01;

f=figure;
for i=1:4;
    FCDvar=FCD+i;
   if i==1
       CD=HCD0;
   elseif i==2
       CD=HCD1;
   elseif i==3
       CD=HCD2;
   else
       CD=HCD3;
   end
[EN,B,G]=CDFQA_energypert_vs_cir_depth(N,HP,HM,CD,psi_initial,Layers,dt,alpha);
 plot(1:Layers,(real(EN(:,1))-D(1,1))/N,'-o','Color',CDcolor{FCDvar},'LineWidth',2)
%plot(1:Layers,(real(EN(:,1))),'-o','Color',CDcolor{FCDvar},'LineWidth',2)

hold on
end

%%
%%


%set(f, 'Units', 'inches', 'Position', [1, 1, 6,4]);  % Set figure size to 6x4 inches

set(f, 'PaperPositionMode', 'auto');
box on
grid on
legend('I','Y','YZ','YX')
%legend('I,Y,YX','YZ','','')
set(gca, "FontName",'Number of Layers','FontSize',24,'LineWidth',1.5)
xlabel('{\cal L}','Interpreter','latex',FontSize=28)
ylabel('e_P',FontSize=28)



exportgraphics(f,'Fig11/GHZ100_pert.jpg')