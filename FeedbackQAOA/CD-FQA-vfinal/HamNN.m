function y=HamNN(N,spin1,spin2,couplingNN,perBC)
 
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

for i=1:N
    if (i==1 || i==N-1)
        q{i}=s1;q{i+1}=s2;
        if i==1
        q(i+2:N)=repmat({speye(2)},1,N-2);
        else
        q(1:N-2)=repmat({speye(2)},1,N-2);   
        end
    elseif (i==N && perBC==1) 
        q{i}=s1;q{1}=s2;q(2:N-1)=repmat({speye(2)},1,N-2);
    else
        q{i}=s1;q{i+1}=s2;
        q(1:i-1)=repmat({speye(2)},1,i-1);
        q(i+2:N)=repmat({speye(2)},1,N-i-1);
    end
    if(i==1)
    H=couplingNN.*mkron(q);
    else
      H=H+couplingNN.*mkron(q);
    end
end
y=H;
end
