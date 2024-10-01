function y=HamNN(N,spin1,spin2,couplingNN,perBC)
usesparse=1;
if usesparse==1

if spin1==0,s11=eye(2);
elseif spin1==1,s11=[0,1;1,0];
elseif spin1==2,s11=[0,-1i;1i,0];
else s11=[1,0;0,-1];
end

if spin2==0,s22=eye(2);
elseif spin2==1,s22=[0,1;1,0];
elseif spin2==2,s22=[0,-1i;1i,0];
else s22=[1,0;0,-1];
end

s1=sparse(s11);
s2=sparse(s22);

H=sparse(zeros(2));
for i=1:N-1
    H=kron(H,sparse(zeros(2)));
end


for i=1:N-1
    if i==1
       
        
            H=H+(couplingNN.*kron(kron(s1,s2),speye(2^(N-2))));
    elseif i==N-1
            H=H+couplingNN.*kron(speye(2^(N-2)),kron(s1,s2));
    else
        H=H+couplingNN.*kron(speye(2^(i-1)),kron(kron(s1,s2),speye(2^(N-2-i+1))));
    end
end
if perBC==1
   H=H+couplingNN.*kron(kron(s2,speye(2^(N-2))),s1);
end
 y=H;
else

%
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
if perBC==1
   H=H+couplingNN.*kron(kron(s2,eye(2^(N-2))),s1);
end
 y=H;
end
end
