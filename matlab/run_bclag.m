clear

x1 = legendre.gauss(100);
x2 = legendre.gauss(2000);
x2(end) = x1(end);



B = bclag_interp_matrix(x1, x2);
B2 = bclag_interp_matrix_mex(x1, x2);

norm(B(:)-B2(:), inf)

timeit(@() bclag_interp_matrix(x1, x2))
timeit(@() bclag_interp_matrix_mex(x1, x2))



x1 = legendre.gauss(10);
x2 = legendre.gauss(20);
x2(1) = x1(1);
x2(end) = x1(end);

B = bclag_interp_matrix(x1, x2);

f = sin(x1/10)
fi = B*f

fi2 = bclag_interp_mex(x1, f, x2)

fi-fi2



