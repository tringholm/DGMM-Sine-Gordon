function out = gradH(u,x)
out = zeros(size(u));
M = length(u)/2;
v = u(M+1:end);
u = u(1:M);
kappa = quadratureWeights(x);
uderiv = 1/4*[(u(2)-u(1))/kappa(1); (u(3:end)-u(1:end-2))./kappa(2:end-1); (u(end)-u(end-1))/kappa(end)];
out(1:M) = kappa.*sin(u);
out(2:M) = out(2:M) + uderiv(1:end-1);
out(1:M-1) = out(1:M-1) - uderiv(2:end);
out(1) = out(1) - uderiv(1);
out(M) = out(M) + uderiv(end);
out(M+1:end) = kappa.*v;
end