clear all

N=6;J=-1;hz=-0.4;hx=-0.4;perBC=1;
CD={'k','b','r','g','m','c','y'};
pool={'I','Y','ZY','YZ'};
FalqonCD=0;% 0 - Original FALQON , 1 - FALOQN with Y, 2 - FALOQN with YZ
Colorcount=0
beta=0;gamma=1;alpha=1;
HP=HamNN(N,3,3,J,perBC)+HamOnsite(N,3,hz)+HamOnsite(N,1,hx);
[EV,D]=eig(HP);

fpopulation0= sprintf('dataCD/FalqonCDpopulation%d.%f.%f.%f.%d.mat',N,J,hz,hx,0);
fpopulation1= sprintf('dataCD/FalqonCDpopulation%d.%f.%f.%f.%d.mat',N,J,hz,hx,1);
fpopulation2= sprintf('dataCD/FalqonCDpopulation%d.%f.%f.%f.%d.mat',N,J,hz,hx,2);
fpopulation3= sprintf('dataCD/FalqonCDpopulation%d.%f.%f.%f.%d.mat',N,J,hz,hx,3);

population0=importdata(fpopulation0);
population1=importdata(fpopulation1);
population2=importdata(fpopulation2);
population3=importdata(fpopulation3);
sz=size(population0);
T=length(population0(1:201));
dt=0.01;


D1=diag(D)-D(1,1);
if mod(ceil(max(D1)),2)==0
    ymax=ceil(max(D1));
else
    ymax=ceil(max(D1))+1;
end

%ymax=ceil(max(D1));
P=zeros((ymax/2)+1,T);
n=1;
     for i=1:length(D1)
         D1(i)
    if (2*(n-1)<= D1(i)) && (D1(i)<2*n)
        P(n,:)=P(n,:)+population3(i,:);
        n
    else
        n=n+1
    end

    % if ((n-1)<= D1(i)) && (D1(i)<n)
    %     P(n,:)=P(n,:)+population3(i,:);
    %     n
    % else
    %     n=n+1
    % end

    end
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

%exportgraphics(f,'3D-FQA3.jpg')

% videofile = VideoWriter('FALQONtfim-012','MPEG-4');
% open(videofile);
% fig1=figure(1);
% ax=gca;
% for i=1:T
%     i
% subplot(1,3,1)
% scatter(diag(D),population0(:,i),100,"filled","o","LineWidth",1)
% hold on
% plot(diag(D),population0(:,i),"LineWidth",4)
% axis([D(1,1) D(2^N,2^N) 0 1])
% str = sprintf('FalqonCD=%s,  dt=%.2f,  Time steps=%d',pool{1},dt,i);
% title(str)
% %title('FalqonCD=',pool{FalqonCD+1},'Time steps=',num2str(i))
% xlabel('Energy levels')
% ylabel('Population')
% fontsize(16,"points")
% time = (i-1)*dt;
% %text(2^(N/2),0.7,num2str(time))
% F=getframe(fig1);
% hold off
% 
% 
% subplot(1,3,2)
% scatter(diag(D),population1(:,i),100,"filled","o","LineWidth",1)
% hold on
% plot(diag(D),population1(:,i),"LineWidth",4)
% axis([D(1,1) D(2^N,2^N) 0 1])
% str = sprintf('FalqonCD=%s,  dt=%.2f,  Time steps=%d',pool{2},dt,i);
% title(str)
% %title('FalqonCD=',pool{FalqonCD+1},'Time steps=',num2str(i))
% xlabel('Energy levels')
% ylabel('Population')
% fontsize(16,"points")
% time = (i-1)*dt;
% %text(2^(N/2),0.7,num2str(time))
% F=getframe(fig1);
% hold off
% 
% 
% subplot(1,3,3)
% scatter(diag(D),population2(:,i),100,"filled","o","LineWidth",1)
% hold on
% plot(diag(D),population2(:,i),"LineWidth",4)
% axis([D(1,1) D(2^N,2^N) 0 1])
% str = sprintf('FalqonCD=%s,  dt=%.2f,  Time steps=%d',pool{3},dt,i);
% title(str)
% %title('FalqonCD=',pool{FalqonCD+1},'Time steps=',num2str(i))
% xlabel('Energy levels')
% ylabel('Population')
% fontsize(16,"points")
% time = (i-1)*dt;
% %text(2^(N/2),0.7,num2str(time))
% F=getframe(fig1);
% hold off
% 
% %pause(0.1)
% % F=getframe(fig1);
% % hold off
% writeVideo(videofile,F);
%end
%poclose(videofile);