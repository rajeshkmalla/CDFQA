function y=HamOnsite(N,spin1,couplingNN)
   
if spin1==0,s1=eye(2);
elseif spin1==1,s1=[0,1;1,0];
elseif spin1==2,s1=[0,-1i;1i,0];
else s1=[1,0;0,-1];
end

H=zeros(2^N);
for i=1:N
    if i==1
           % H=H-hz.*blkdiag(sz,eye(2^(N-1)));
            H=H+kron(s1,eye(2^(N-1)));
    elseif i==N
            %H=H-hz.*blkdiag(eye(2^(N-1)),sz);
            H=H+kron(eye(2^(N-1)),s1);
    else
       % H=H-hz.*blkdiag(eye(2^(i-1)),sz,eye(2^(N-i)));
        H=H+kron(eye(2^(i-1)),kron(s1,eye(2^(N-i))));
    end
end

 y=couplingNN.*H;
end
