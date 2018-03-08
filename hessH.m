function out = hessH(u_new,x)
M = length(u_new)/2;
u_new = u_new(1:M);
O = sparse(M,M);
kappa = quadratureWeights(x);
i = 1:M;
B = sparse(i,i,kappa,M,M);
% B = spdiags(kappa,0,M,M);

Asupdiag = -1./(4*kappa(2:end-1));
Adiag = 1/4*(1./([kappa(1); kappa(1:end-1)]) + 1./([kappa(2:end); kappa(end)])) + kappa.*cos(u_new);
Asubdiag = -1./(4*kappa(2:end-1));

i = [3:M, 1:M, 1:M-2];
j = [1:M-2, 1:M, 3:M];
A = sparse(i,j,[Asubdiag; Adiag; Asupdiag],M,M);

% A = spdiags([Asubdiag Adiag Asupdiag],-2:2:2,M,M);
A(1,2) = -1/(4*kappa(1));
A(2,1) = -1/(4*kappa(1));
A(end,end-1) = -1/(4*kappa(end));
A(end-1,end) = -1/(4*kappa(end));
out = [A O;
       O B];
end