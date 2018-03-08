function kappa = quadratureWeights(x)
kappa = zeros(size(x));
kappa(1) = 1/2*(x(2)-x(1));
kappa(2:end-1) = 1/2*(x(3:end)-x(1:end-2));
kappa(end) = 1/2*(x(end)-x(end-1));
end