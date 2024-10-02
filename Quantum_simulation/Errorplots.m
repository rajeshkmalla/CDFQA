dataFCDy = [-1.8392700229319439, -2.2666712462333316, ...
-2.936265422867532, -3.732462126362557, -4.335154531302772, ...
-4.59504045908713, -4.663320886131266, -4.676413580330181, ...
-4.678384780480009, -4.6785644797493005];
dataFCD300 = [-1.92533333, -2.25066667, -2.81333333, -3.61866667, ...
-4.45866667, -4.91066667, -4.71066667, -4.63466667, -4.67866667, ...
-4.65466667];
dataFCD300std = [0.10998434905493162, 0.10999318497407345, ...
   0.10787400057474461, 0.09994078246641859, 0.0809776558115925, ...
   0.06164687970836624, 0.061302165993996946, 0.060433691854056976, ...
   0.058323796045740936, 0.06382415657854762];
dataFCD3000 = [-1.7752, -2.23693333, -2.9272, -3.78, -4.29946667, ...
-4.57386667, -4.68466667, -4.63133333, -4.67386667, -4.7016];
dataFCD3000std = [0.034821765011534354, 0.03477581151701088, ...
   0.03395317839264777, 0.0308792859253458, 0.026703185439055633, ...
   0.022476278590088248, 0.020142792627012143, 0.019755182648648974, ...
   0.019039390902177673, 0.018914680639749215];
dataFCD30000 = [-1.85654667, -2.27029333, -2.90276, -3.73830667, ...
-4.32369333, -4.59130667, -4.66542667, -4.69076, -4.66926667, ...
-4.68890667];
dataFCD30000std = [0.011014023374401011, 0.010986639077095888, ...
   0.010755295854807737, 0.009823111524849427, 0.0083601031824219, ...
   0.007082278674111729, 0.0063924062852384125, 0.006091747445945637,... 
   0.006054848499068054, 0.005985016717944411];
dataFCD300000 = [-1.837912, -2.27214133, -2.94936667, -3.74052133, ...
-4.33484533, -4.59592, -4.66627733, -4.68266267, -4.676156, ...
-4.67461333];
dataFCD300000std = [0.003483767699706316, 0.0034739982565594038,... 
   0.0033900260954948854, 0.003106598235424675, 0.002634154254475872, ...
   0.0022307111991012097, 0.0020239986539125963, ...
   0.0019357055542158841, 0.0019134150846476146, ...
   0.001903216023151145];
N=4;
f=figure;
E0=-4.78086132;
% plot(1:length(dataFCDy),(dataFCDy-E0)./N,'Color','k','LineWidth',2)
% hold on
% e1=errorbar(1:length(dataFCD300),(dataFCD300-E0)./N,dataFCD300std./N,'LineWidth',2)
% e2=errorbar(1:length(dataFCD3000),(dataFCD3000-E0)./N,dataFCD3000std./N,'LineWidth',2)
% e3=errorbar(1:length(dataFCD30000),(dataFCD30000-E0)./N,dataFCD30000std./N,'LineWidth',2)
% e4=errorbar(1:length(dataFCD300000),(dataFCD300000-E0)./N,dataFCD300000std./N,'LineWidth',2)
% e1.Color = '#0072BD';
%  e2.Color='#D95319';
%  e3.Color='#7E2F8E';
%  e4.Color='m';

% 
inv_samplesize=1./sqrt([300,3000,30000,300000]);
x=inv_samplesize;
x0=1./sqrt([300,3000,30000,300000,30000000000000]);
variance_data=[0.393898, 0.109965, 0.0448779, 0.0183739]./N;
y1=variance_data;
P = polyfit(inv_samplesize,variance_data,1);
yfit = polyval(P,x0);
scatter(inv_samplesize,variance_data,100,"filled")
hold on
plot(x0,yfit,'LineWidth',2);
box on
legend('','y=1.69 x')
% eqn = string(" Linear: y = " + "1.69") + "x  ";
% text(min(x),max(y1),eqn,"HorizontalAlignment","left","VerticalAlignment","top",FontSize=18)
%  hold on
% e11=errorbar(1:length(dataFCD300),(dataFCD300-dataFCDy)./N,dataFCD300std./N,'LineWidth',2)
% e21=errorbar(1:length(dataFCD3000),(dataFCD3000-dataFCDy)./N,dataFCD3000std./N,'LineWidth',2)
% e31=errorbar(1:length(dataFCD30000),(dataFCD30000-dataFCDy)./N,dataFCD30000std./N,'LineWidth',2)
% e41=errorbar(1:length(dataFCD300000),(dataFCD300000-dataFCDy)./N,dataFCD300000std./N,'LineWidth',2)
% e11.Color = '#0072BD';
% e21.Color='#D95319';
% e31.Color='#7E2F8E';
% e41.Color='m';
box on
 grid on
 
%legend('C','300','3000','30000','300000')
set(gca, "FontName",'Number of Layers','FontSize',24,'LineWidth',1.5)
% xlabel('$1/{\sqrt M}$','Interpreter','latex',FontSize=28)
% ylabel('(\Delta e_P)^2',FontSize=28)
% 1/\sqrt{M}
xlabel('$1/\sqrt{M}$','Interpreter','latex',FontSize=28)
ylabel('(\Deltae_P)^2',FontSize=28)
%axis([1 10 -0.1 0.1])

% e11=errorbar(1:length(dataFCD300),(dataFCD300-dataFCDy)./N,dataFCD300std./N,'LineWidth',2)
% e21=errorbar(1:length(dataFCD3000),(dataFCD3000-dataFCDy)./N,dataFCD3000std./N,'LineWidth',2)
% e31=errorbar(1:length(dataFCD30000),(dataFCD30000-dataFCDy)./N,dataFCD30000std./N,'LineWidth',2)
% e41=errorbar(1:length(dataFCD300000),(dataFCD300000-dataFCDy)./N,dataFCD300000std./N,'LineWidth',2)
% e11.Color = '#0072BD';
% e21.Color='#D95319';
% e31.Color='#7E2F8E';
% e41.Color='m';
