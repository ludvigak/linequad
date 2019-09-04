function A = full_matrix(g)

N = numel(g.z);
A = zeros(N,N);
disp(' * Assembly and solve')
tic
for i=1:N
    for j=1:N
        A(i,j) = imag( g.dz(j)/(g.z(j)-g.z(i)) )/(2*pi);
    end
    A(i,i) = 1/2 + g.w(i).*imag(g.zpp(i)./(2*g.zp(i)))/(2*pi);
end
