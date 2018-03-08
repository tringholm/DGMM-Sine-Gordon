function out = F_DVDdir(u_new,u_int,diffH,dt,x)
M = length(u_new)/2;
I = speye(M,M);
O = sparse(M,M);
S = [O I;
    -I O];
DVD = DVD_H(u_new,u_int,x);

kappa = quadratureWeights(x);
kappa2 = [kappa; kappa];

phi = (kappa2.*DVD)'*DVD;
out = 1/dt*(u_new-u_int) - diffH/(phi*dt)*DVD - S*DVD;

end