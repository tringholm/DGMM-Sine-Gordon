function out = F_gradHdir(u_new,u_int,diffH,dt,x)
M = length(u_new)/2;
I = speye(M,M);
O = sparse(M,M);
S = [O I;
    -I O];
grad = gradH(u_new,x);
DVD = DVD_H(u_new,u_int,x);

kappa = quadratureWeights(x);
kappa2 = [kappa; kappa];

phi = (kappa2.*DVD)'*grad;
out = 1/dt*(u_new-u_int) - diffH/dt*grad/phi - S*DVD;
end