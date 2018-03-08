function out = jac_DVD_H(u_new,u_int,x)
kappa = quadratureWeights(x);
M = length(u_new)/2;
u_new = u_new(1:M);
u_int = u_int(1:M);
I  = speye(M);
O = sparse(M,M);

bad_indices = ~(u_new-u_int);
good_indices = ~bad_indices;
Adiag_trigterms = zeros(size(u_new));
Adiag_trigterms(good_indices) = sin(u_new(good_indices))./(u_new(good_indices)-u_int(good_indices))...
                             -2*sin((u_new(good_indices)+u_int(good_indices))/2).*sin((u_new(good_indices)-u_int(good_indices))/2)./((u_new(good_indices)-u_int(good_indices)).^2);
Adiag_trigterms(bad_indices) = 1/2*cos(u_int(bad_indices));
Asupdiag = -1/8*1./(kappa(1:end-2).*kappa(2:end-1));
Adiag = 1./(8*kappa).*(1./[kappa(1); kappa(1:end-1)] + 1./([kappa(2:end); kappa(end)])) + Adiag_trigterms;
Asubdiag = -1/8*1./(kappa(3:end).*kappa(2:end-1));

i = [3:M, 1:M, 1:M-2];
j = [1:M-2, 1:M, 3:M];
A = sparse(i,j,[Asubdiag; Adiag; Asupdiag],M,M);

% A = spdiags([Asubdiag Adiag Asupdiag],-1:1,M,M);
A(1,2) = -1/(8*kappa(1)^2);
A(2,1) = -1/(8*kappa(1)*kappa(2));
A(end-1,end) = -1/(8*kappa(end-1)*kappa(end));
A(end,end-1) = -1/(8*kappa(end)*kappa(end));
out = [A O;
       O 1/2*I];
end