function y=HamTVCD(N,couplingNN)
s1=[0,1;1,0];
s2=[0,-1i;1i,0];
s3=[1,0;0,-1];
H=zeros(4^N);

for i=1:N-1       
        if i==1
            H=H+couplingNN.*kron(kron((kron(s1,s2)-kron(s2,s1)),eye(2^(N-2))),kron(s3,eye(2^(N-1)))-kron(eye(2),kron(s3,eye(2^(N-2)))));
    elseif i==N-1
            H=H+couplingNN.*kron(kron(eye(2^(N-2)),(kron(s1,s2)-kron(s2,s1))),kron(kron(eye(2^(N-2)),s3),eye(2))-kron(eye(2^(N-1)),s3));
    else
        H=H+couplingNN.*kron(kron(eye(2^(i-1)),kron((kron(s1,s2)-kron(s2,s1)),eye(2^(N-2-i+1)))),kron(eye(2^(i-1)),kron(s3,eye(2^(N-i))))-kron(eye(2^(i)),kron(s3,eye(2^(N-i-1)))));
        end
end

 y=0.25*H;
end