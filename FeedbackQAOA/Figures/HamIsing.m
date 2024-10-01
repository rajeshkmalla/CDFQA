function y=HamIsing(N,J,hz,hx)
H=zeros(2^N);

%%%%%%%%%%%%%
I=eye(2);
sz=[1,0;0,-1];
sx=[0,1;1,0];
sy=[0,-1i;1i,0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:N-1
    if i==1
            H=H-J.*kron(kron(sz,sz),eye(2^(N-2)));
    elseif i==N-1
            H=H-J.*kron(eye(2^(N-2)),kron(sz,sz));
    else
        H=H-J.*kron(eye(2^(i-1)),kron(kron(sz,sz),eye(2^(N-2-i+1))));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:N
    if i==1
            H=H-hz.*kron(sz,eye(2^(N-1)));
            H=H-hx.*kron(sx,eye(2^(N-1)));
    elseif i==N
            H=H-hz.*kron(eye(2^(N-1)),sz);
            H=H-hx.*kron(eye(2^(N-1)),sx);
    else
        H=H-hz.*kron(eye(2^(i-1)),kron(sz,eye(2^(N-i))));
        H=H-hx.*kron(eye(2^(i-1)),kron(sx,eye(2^(N-i))));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y=H;
end