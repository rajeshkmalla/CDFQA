clear all
close all
%% Initialization
N=6;J=-1;hz=-0.4;hx=-0;perBC=1;% Hamiltonian parameters
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
Layers=800;
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
[EN,B,G]=CDFQA_energy_vs_cir_depth(N,HP,HM,CD,psi_initial,Layers,dt,alpha);
%plot(1:Layers,(real(EN(:,1))-D(1,1))/N,'-o','Color',CDcolor{FCDvar},'LineWidth',2)
%plot(1:Layers,(real(EN(:,1))),'-o','Color',CDcolor{FCDvar},'LineWidth',2)

hold on

mdata(:,FCDvar)=(real(EN(:,1))-D(1,1))/N;

for t=1:Layers

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

%%


    plot(1:Layers/2,Xdata(1:Layers/2,1),'o','Color',CDcolor{1},'LineWidth',2)
    hold on
    plot(1:Layers,Xdata(1:Layers,2),'o','Color',CDcolor{2},'LineWidth',2)
    plot(1:Layers,Xdata(1:Layers,3),'o','Color',CDcolor{3},'LineWidth',2)
    plot(1:Layers*2,Xdata(1:Layers*2,4),'o','Color',CDcolor{4},'LineWidth',2)
    
    h1=plot(Layers/2:Layers*2,Xdata(Layers/2:Layers*2,1),'o','Color',CDcolor{1},'LineWidth',2)
    hold on
    h2=plot(Layers:Layers*2,Xdata(Layers:Layers*2,2),'o','Color',CDcolor{2},'LineWidth',2)
    h3=plot(Layers:Layers*2,Xdata(Layers:Layers*2,3),'o','Color',CDcolor{3},'LineWidth',2)
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
axis([0.01 Layers*2 0.01 1.5])


%%
%%

% 
% %set(f, 'Units', 'inches', 'Position', [1, 1, 6,4]);  % Set figure size to 6x4 inches
% 
% set(f, 'PaperPositionMode', 'auto');
% box on
% grid on
% %legend('I','Y','YZ','YX')
% %legend('I,Y','YZ')
% set(gca, "FontName",'Number of Layers','FontSize',24,'LineWidth',1.5)
% xlabel('{\cal L}','Interpreter','latex',FontSize=28)
% ylabel('e_P',FontSize=28)
% 
% %print(fig, 'figure_output', '-dpng', '-r300');  % PNG file with 300 dpi resolution
 exportgraphics(f,'Fig3/measurements.jpg')