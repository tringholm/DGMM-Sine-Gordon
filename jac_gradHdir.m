function out = jac_gradHdir(u_new,u_int,diffH,dt,x)
M = length(u_new)/2;
I = speye(M,M);
O = sparse(M,M);
S = [O I;
    -I O];

kappa = quadratureWeights(x);
kappa2 = [kappa; kappa];

grad = gradH(u_new,x);
DVD = DVD_H(u_new,u_int,x);
hess = hessH(u_new,x);
DVDjac = jac_DVD_H(u_new,u_int,x);

phi = (kappa2.*DVD)'*grad;
gradphi = DVDjac'*(kappa2.*grad) + hess*(kappa2.*DVD);

out = 1/dt*speye(2*M,2*M) - S*DVDjac - diffH/dt*(1/phi*hess - 1/phi^2*grad*gradphi');
end