function y=HM(N)
H=zeros(2^N);
I=eye(2);
sz=[1,0;0,-1];
sx=[0,1;1,0];
sy=[0,-1i;1i,0];
for i=1:N
    if i==1
           % H=H-hz.*blkdiag(sz,eye(2^(N-1)));
            H=H+kron(sx,eye(2^(N-1)));
    elseif i==N
            %H=H-hz.*blkdiag(eye(2^(N-1)),sz);
            H=H+kron(eye(2^(N-1)),sx);
    else
       % H=H-hz.*blkdiag(eye(2^(i-1)),sz,eye(2^(N-i)));
        H=H+kron(eye(2^(i-1)),kron(sx,eye(2^(N-i))));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y=H;
end