function y=Creation(m,N)
s1=[0,1;1,0];
s2=[0,-1i;1i,0];
s3=[1,0;0,-1];

H=zeros(2^N);
for i=1:N
    
    if i==1
           % H=H-hz.*blkdiag(sz,eye(2^(N-1)));
            H=H+kron(s1-1i*s2,eye(2^(N-1)));
            s3z=s3;
    elseif i==N
            %H=H-hz.*blkdiag(eye(2^(N-1)),sz);
            H=H+kron(s3z,s1-1i*s2);
    else
       % H=H-hz.*blkdiag(eye(2^(i-1)),sz,eye(2^(N-i)));
        H=H+kron(s3z,kron(s1-1i*s2,eye(2^(N-i))));
         s3z=kron(s3z,s3);
    end
end

 y=H;
end
