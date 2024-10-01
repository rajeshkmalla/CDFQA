function y=HamNN(N,spin1,spin2,couplingNN)

 
%getopt('INIT',varargin);
% spin1=getopt('spin1',1);
% spin2=getopt('spin2',1);
% couplingNN=getopt('couplingNN',1);
   
if spin1==0,s1=eye(2);
elseif spin1==1,s1=[0,1;1,0];
elseif spin1==2,s1=[0,-1i;1i,0];
else s1=[1,0;0,-1];
end
    
if spin2==0,s2=eye(2);
elseif spin2==1,s2=[0,1;1,0];
elseif spin2==2,s2=[0,-1i;1i,0];
else s2=[1,0;0,-1];
end

H=zeros(2^N);
for i=1:N-1
    if i==1
            H=H+couplingNN.*kron(kron(s1,s2),eye(2^(N-2)));
    elseif i==N-1
            H=H+couplingNN.*kron(eye(2^(N-2)),kron(s1,s2));
    else
        H=H+couplingNN.*kron(eye(2^(i-1)),kron(kron(s1,s2),eye(2^(N-2-i+1))));
    end
end
 y=H;
end
