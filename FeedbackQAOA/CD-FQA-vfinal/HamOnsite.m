function y=HamOnsite(N,spin1,couplingonsite)
   

if spin1==0,s11=eye(2);
elseif spin1==1,s11=[0,1;1,0];
elseif spin1==2,s11=[0,-1i;1i,0];
else s11=[1,0;0,-1];
end

s1=sparse(s11);

for i=1:N
    if i==1
        q{i}=s1;
        q(i+1:N)=repmat({speye(2)},1,N-1);    
    elseif i==N
        q{i}=s1;
        q(1:N-1)=repmat({speye(2)},1,N-1);  
    else
       q{i}=s1;
       q(1:i-1)=repmat({speye(2)},1,i-1);
       q(i+1:N)=repmat({speye(2)},1,N-i);
    end
    
    if(i==1)
    H=mkron(q);
    else
      H=H+mkron(q);
    end
end
 y=couplingonsite.*H;
end
