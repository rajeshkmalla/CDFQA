function y=mkron(q)
y=q{1};
for i=2:length(q),y=kron(y,q{i});end
end