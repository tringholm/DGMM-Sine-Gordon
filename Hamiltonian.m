function out = Hamiltonian(u,x)
n = length(u)/2;
v = u(n+1:end);
u = u(1:n);
kappa = quadratureWeights(x);
u_derivs = zeros(length(u),1);
u_derivs(1) = (u(2)-u(1))/kappa(1);
u_derivs(2:end-1) = (u(3:end)-u(1:end-2))./kappa(2:end-1);
u_derivs(end) = (u(end)-u(end-1))/kappa(end);
sum1 = 1/2*sum(sort(kappa.*v.^2));
sum2 = 1/8*sum(sort(kappa.*u_derivs.^2));
sum3 = -sum(sort(kappa.*cos(u)));
sum4 = x(end)-x(1);
out = sort([sum1 sum2 sum3 sum4]);
out = sum(out);
end