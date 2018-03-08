function [f, J] = F(u,u_int,diffH,dt,x_new,projtype)
% Function for evaluating system of equations to be solved at each time
% step, and the Jacobian of this.
if projtype == 1 % Correct for energy difference in direction of gradH 
    J = jac_gradHdir(u,u_int,diffH,dt,x_new);
    f = F_gradHdir(u,u_int,diffH,dt,x_new);
else % Correct for energy difference in direction of DVD 
    J = jac_DVDdir(u,u_int,diffH,dt,x_new);
    f = F_DVDdir(u,u_int,diffH,dt,x_new);
end
end