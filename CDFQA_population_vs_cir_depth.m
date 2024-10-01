function population=cdfqa_population_vs_cir_depth(N,HP,HM,HCD,psi_initial,Layers,dt,alpha)

%%
%EN=zeros(Layers,1);  
B=zeros(Layers,1);  
G=zeros(Layers,1); 
population=zeros(2^N,Layers);

CM=(HM*HP-HP*HM);
CMCD=(HCD*(HP))-((HP)*HCD);
[EV,D]=eig(HP);
%%
B(1)=0;G(1)=0;EN(1)=(psi_initial'*HP*psi_initial);psi=psi_initial;
population(:,1)=abs((psi'*EV)').^2;
%%
for i=1:Layers-1
    t=dt.*(i-1);
    %
    psi=fastExpm(-1i*dt*G(i)*HCD)*(fastExpm(-1i*dt*B(i)*HM)*(fastExpm(-1i*dt*HP)*psi));
    %
    B(i+1)=-1i.*alpha*(psi'*CM*psi);
    G(i+1)=-1i.*alpha*(psi'*CMCD*psi);
    %EN(i+1)=(psi'*HP*psi);
    population(:,i+1)=abs((psi'*EV)').^2;
    %
end
end
