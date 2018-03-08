function out = DVD_H(u_new,u_int,x)
M=length(u_new)/2;
out = zeros(size(u_new));
v_new = u_new(M+1:end);
v_int = u_int(M+1:end);
u_new = u_new(1:M);
u_int = u_int(1:M);
kappa = quadratureWeights(x);
unewderiv = 1/8*[(u_new(2)-u_new(1))/kappa(1); (u_new(3:end)-u_new(1:end-2))./kappa(2:end-1); (u_new(end)-u_new(end-1))/kappa(end)];
uintderiv = 1/8*[(u_int(2)-u_int(1))/kappa(1); (u_int(3:end)-u_int(1:end-2))./kappa(2:end-1); (u_int(end)-u_int(end-1))/kappa(end)];

bad_indices = ~(u_new-u_int);
good_indices = ~bad_indices;

out(good_indices) = 2*sin((u_new(good_indices)+u_int(good_indices))/2).*sin((u_new(good_indices)-u_int(good_indices))/2)./(u_new(good_indices)-u_int(good_indices));
out(bad_indices) = sin(u_int(bad_indices));
out(2:M) = out(2:M) + 1./(kappa(2:end)).*(unewderiv(1:end-1) + uintderiv(1:end-1));
out(1:M-1) = out(1:M-1) - 1./(kappa(1:end-1)).*(unewderiv(2:end) + uintderiv(2:end));
out(1) = out(1) - (unewderiv(1) + uintderiv(1))/kappa(1);
out(end) = out(end) + (unewderiv(end) + uintderiv(end))/kappa(end);
out(M+1:end) = 1/2*(v_new+v_int);
end