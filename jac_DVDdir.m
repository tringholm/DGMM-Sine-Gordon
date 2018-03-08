function out = jac_DVDdir(u_new,u_int,diffH,dt,x)
M = length(u_new)/2;
I = speye(M,M);
O = sparse(M,M);
S = [O I;
    -I O];

kappa = quadratureWeights(x);
kappa2 = [kappa; kappa];

DVD = DVD_H(u_new,u_int,x);
DVDjac = jac_DVD_H(u_new,u_int,x);

phi = (kappa2.*DVD)'*DVD;


out = 1/dt*speye(2*M,2*M) - S*DVDjac - diffH/dt*(1/phi*DVDjac - 2/phi^2*DVD*(kappa2.*DVD)'*DVDjac);

end