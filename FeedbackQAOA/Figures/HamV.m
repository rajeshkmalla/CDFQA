function y=HamV(N,couplingU)
s3=[1,0;0,-1];
H=zeros(4^N);





for i=1:N  
    if i==1
            H=H+couplingU.*kron(kron((eye(2)-s3),eye(2^(N-1))),kron((eye(2)-s3),eye(2^(N-1))));
    elseif i==N
           H=H+couplingU.*kron(kron(eye(2^(N-1)),(eye(2)-s3)),kron(eye(2^(N-1)),(eye(2)-s3)));
    else
            H=H+couplingU.*kron(kron(eye(2^(i-1)),kron((eye(2)-s3),eye(2^(N-i)))),kron(eye(2^(i-1)),kron((eye(2)-s3),eye(2^(N-i)))));
    end
end
 y=-0.25.*H;
end