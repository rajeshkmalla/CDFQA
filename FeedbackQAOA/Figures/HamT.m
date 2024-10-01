function y=HamT(N,couplingNN)
s1=[0,1;1,0];
s2=[0,-1i;1i,0];
H=zeros(2^N);

for i=1:N-1        
        if i==1
            H=H+couplingNN.*kron((kron(s1,s1)+kron(s2,s2)),eye(2^(N-2)));
    elseif i==N-1
            H=H+couplingNN.*kron(eye(2^(N-2)),(kron(s1,s1)+kron(s2,s2)));
    else
        H=H+couplingNN.*kron(eye(2^(i-1)),kron((kron(s1,s1)+kron(s2,s2)),eye(2^(N-2-i+1))));
        end
end
 y=-0.5.*kron(H,H);
end
