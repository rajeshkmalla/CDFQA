function [EN,B,G]=CDFQA_energy_vs_cir_depth(N,HP,HM,HCD,psi_initial,Layers,dt,alpha)

%%
EN=zeros(Layers,1);  
B=zeros(Layers,1);  
G=zeros(Layers,1); 
%population=zeros(2^N,2^N);

CM=(HM*HP-HP*HM);
CMCD=(HCD*(HP))-((HP)*HCD);
%%
B(1)=0;G(1)=0;EN(1)=(psi_initial'*HP*psi_initial);psi=psi_initial;
%population(:,1)=abs(psi'*EV)').^2;
%%
for i=1:Layers
    t=dt.*(i-1);
    %
    B(i)=-1i.*alpha*(psi'*CM*psi);
    G(i)=-1i.*alpha*(psi'*CMCD*psi);
    EN(i)=(psi'*HP*psi);
    psi=fastExpm(-1i*dt*G(i)*HCD)*(fastExpm(-1i*dt*B(i)*HM)*(fastExpm(-1i*dt*HP)*psi));
    %population(:,i+1)=abs((psi'*EV)').^2;
    %
end
end
