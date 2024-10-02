clear all
close all
%% Initialization
N=6;J=-1;hz=-0.4;hx=-0.4;perBC=1;% Hamiltonian parameters
CD={"#7E2F8E",'k',"#0072BD","#D95319",'g','c','y'};
pool={'I','Y','ZY','YZ'};
FalqonCD=0;% 0-original FALQON , 1 - FALOQN with Y, 2 - FALOQN with YZ

BETA=[2,0.1,0.5,1];%TIME=[0.005,0.01,0.02,0.04]./2;
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

f=figure;
Layers=1000;
dt=0.01;
for k=1:4
    alpha=BETA(k)*(6/N);
    %k1=dt/0.02;
    %Layers=10/dt;
   [EN,B,G]=CDFQA_energy_vs_cir_depth(N,HP,HM,HCD3,psi_initial,Layers,dt,alpha);
loglog((1:Layers),(real(EN(:,1))-D(1,1))/N,'Color',CD{k},'LineWidth',2)
hold on
grid on
end


%set(f, 'Units', 'inches', 'Position', [1, 1, 6,4]);  % Set figure size to 6x4 inches

%set(f, 'PaperPositionMode', 'auto');
box on
grid on
%legend('I','Y','YZ','YX')
%legend('I,Y','YZ')
set(gca, "FontName",'Number of Layers','FontSize',24,'LineWidth',1.5)
% xlabel('{\cal L}','Interpreter','latex',FontSize=28)
% ylabel('e_P',FontSize=28)

%print(fig, 'figure_output', '-dpng', '-r300');  % PNG file with 300 dpi resolution
exportgraphics(f,'Fig6/MFIdifferentalpha3.jpg')