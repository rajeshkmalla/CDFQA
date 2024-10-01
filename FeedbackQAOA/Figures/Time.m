function y=Time(N,J,hz,hx)
dt = 0.01;
tmin = 0; 
tmax = 10;
L=length(tmin:dt:tmax);
v=zeros(2^N,L);
v(:,1)=ones(2^N,1)./sqrt(2^N);
for i=1:L-1
    t=tmin+(dt*(i-1));
    h=Ham;
   v(:,i+1)=expm(-1i*dt*h)*v(:,i);
end
y=v;
end