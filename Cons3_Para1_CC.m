r1 = (alpha1r*C01(1)/(k+1));
r2 = (alpha2r*C01(1)/(k+1));
d1 = -alpha1d*C01(2);
d2 = -alpha2d*C01(2);
E_cons = 0.2762/2 + (((k+1)^2)*r1*r2 + d1*d2 + (k+1)*r1*d1 + (k+1)*r2*d2 - (k+1)*r2*d1 - (k+1)*r1*d2);
Ac = ((k+1)^2)*r1*r2 - (k+1)*r1*d2;
Bc = (d1*d2 - (k+1)*r2*d1);
Cc = (k+1)*r1*d1 + (k+1)*r2*d2;
Ep1 = E_cons + 2*Cc;
deltest = pi/4;
delval = 4*Bc*(cos(deltest)^3) + 2*Cc*(cos(deltest)^2) + (Ac - 3*Bc)*(cos(deltest));

checkmat = [Ep1 delval Ep1-delval]