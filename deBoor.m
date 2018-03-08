function x_new = deBoor(x,u)
kappa = 10;
M = length(x)-1;
xdiff = x(2:end)-x(1:end-1);
p = sqrt(1 + kappa*((u(2:end) - u(1:end-1))./xdiff).^2);
p(1) = 0.5*(p(1)+p(2));
p(end) = 0.5*(p(end-1) + p(end));
p(2:end-1) = 0.25*p(1:end-2) + 0.5*p(2:end-1) + 0.25*p(3:end);
P = [0; cumsum(p.*xdiff)];
p = [0; p];
x_new = zeros(size(x));
x_new(1) = x(1);
x_new(end) = x(end);
for i = 2:M
   disc = (i-1)/M*P(end);
   k = find(~(disc > P),1);
   x_new(i) = x(k-1) + ((i-1)/M*P(end) - P(k-1))/p(k); 
end
end