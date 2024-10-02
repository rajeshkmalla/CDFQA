count=0;
f=figure;
depth=zeros(9,1);
for N=4:12
    count=count+1;
fenergy= sprintf('Plots/N.%d.mat',N);
energylog=real(importdata(fenergy));
%energylog=log10(energy);
sz=length(energylog);
for i=1:sz
    if(floor(energylog(i))==-1 && floor(energylog(i+1))==-2)
    depth(count,1)=i+1
    end
end
end
loglog(4:12,depth,'Color','#0072BD','LineWidth',3)
hold on
grid on
%scatter(4:12,depth,100,"filled",'Color','red')
%legend(CD1{:})
set(gca, "FontName",'Number of Layers','FontSize',24,'LineWidth',1.5)
xlabel('{\cal N}','Interpreter','latex',FontSize=28)
ylabel('{\cal L}','Interpreter','latex',FontSize=28)
exportgraphics(f,'Plots/errorplot.jpg')