%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IF USING THIS CODE FOR RESEARCH PURPOSES, PLEASE CITE OUR ARTICLE     %
% Eidnes, S., Owren, B. & Ringholm, T. Adv Comput Math (2017).          %
% https://doi.org/10.1007/s10444-017-9562-8                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

nx = 401;       % Number of spatial discretization points
nt = 201;       % Number of time steps
xmin = -30;     % Left edge of interval
xmax = 30;      % Right edge of interval
tmin = 0;       % Starting time
tmax = 10;       % Stopping time
moving = 1;     % True if using moving mesh
doplot = 1;     % True if producing plots while running
projtype = 2;   % Projection direction: 1 for grad H, otherwise DVD direction
interptype = 2; % Interpolation type: 1 for linear, otherwise cubic

eqntype = 2;    % Three demo equations

if eqntype == 1     % Breather solution
    omega = 0.2;
    u_analytic = @(x,t) 4*atan((sqrt(1-omega^2)*sin(omega*t))./(omega*cosh(sqrt(1-omega^2)*x)));
    ux_analytic = @(x,t) -sqrt(1-omega^2)*tanh(sqrt(1-omega^2)*x).*u_analytic(x,t);
    v_analytic = @(x,t) 4*sqrt(1-omega^2)*cos(omega*t)./cosh(sqrt(1-omega^2)*x).*1./(1+(1-omega^2)*sin(omega*t)^2./(omega^2*(cosh(sqrt(1-omega^2)*x).^2)));
elseif eqntype == 2     % Kink-antikink solution
    epsilon = 0.01;
    c = 1-epsilon;
    u_analytic = @(x,t) 4*atan(sinh(c*t/sqrt(1-c^2))./(c*cosh(1/sqrt(1-c^2)*x)));
    v_analytic = @(x,t) 4*c^2/sqrt(1-c^2)*cosh(c*t/sqrt(1-c^2))*cosh(1/sqrt(1-c^2)*x)./(c^2*cosh(1/sqrt(1-c^2)*x).^2 + sinh(c*t/sqrt(1-c^2))^2);
else     % Kink solution
    c = 0.8;
    x_0 = -10;
    u_analytic = @(x,t) 4*atan(exp(1/sqrt(1-c^2)*(x-x_0-c*t)));
    v_analytic = @(x,t) -2*c/sqrt(1-c^2)*1./(cosh(1/sqrt(1-c^2)*(x-x_0-c*t)));
end

% Initial grid and function values
dx = (xmax-xmin)/(nx-1);
dt = (tmax-tmin)/(nt-1);
x = xmin:dx:xmax;
x = x';
u0 = @(z) u_analytic(z,tmin);
v0 = @(z) v_analytic(z,tmin);

% Run function with plotting of analytic solutions
[u_collection, x_collection, H] = SineGordonAVF(u0,v0,x,dt,tmin,tmax,moving,projtype,interptype,doplot,u_analytic,v_analytic);

% Run function without plotting of analytic solutions
[u_collection, x_collection, H] = SineGordonAVF(u0,v0,x,dt,tmin,tmax,moving,projtype,interptype,doplot);

% Run function without moving mesh or plotting
moving = 0;
doplot = 0;
[u_collection, x_collection, H] = SineGordonAVF(u0,v0,x,dt,tmin,tmax,moving,projtype,interptype,doplot);



