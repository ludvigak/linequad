function new = diff_curve(curve)

d = diff(curve);
d2 = diff(d);

new.sym = curve;
new.tau = matlabFunction(curve);
new.dtau = matlabFunction(d);
new.d2tau = matlabFunction(d2);
new.normal = matlabFunction(1i*d / abs(d));
new.tangent = matlabFunction(d / abs(d));
