function [u_collection,x_collection,H] = SineGordonAVF(u0_func,v0_func,x_old,dt,tmin,tmax,moving,projtype,interptype,doplot,u_analytic,v_analytic)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IF USING THIS CODE FORE RESEARCH PURPOSES, PLEASE CITE OUR ARTICLE    %
% Eidnes, S., Owren, B. & Ringholm, T. Adv Comput Math (2017).          %
% https://doi.org/10.1007/s10444-017-9562-8                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve the 1D Sine-Gordon equation numerically using moving mesh and the
% AVF scheme for time stepping.
%
% Input: 
%
% u0            - inline function for evaluating initial condition on u
% v0            - inline funciton for evaluating initial condition on v = u'
% x             - initial mesh
% dt            - time step size
% tmin          - starting time 
% tmax          - stopping time
% moving        - set to true to use moving mesh
% projtype      - choose correction direction for moving mesh
% interptype    - choose interpolation type. 1 means linear, otherwise cubic
% doplot        - set to true for plotting while function runs
% u_analytic    - analytic solution for u, only for testing purposes
% v_analytic    - analytic solution for v, only for testing purposes
%
% Output: 
%
% u_collection  - collection of solution u at all time steps
% x_collection  - collection of all meshes used
% H             - Hamiltonian values at all time steps
%



if nargin > 11
    have_u_analytic = 1;
else
    have_u_analytic = 0;
end
if nargin > 12
    have_v_analytic = 1;
else
    have_v_analytic = 0;
end
x_new = x_old;

N = length(tmin:dt:tmax);
M = length(x_old);

u0 = u0_func(x_old);
v0 = v0_func(x_old);

% Initialize grid
change = inf;
while change > 0.001 && moving
    x_new = deBoor(x_old,u0);
    u0 = u0_func(x_new);
    v0 = v0_func(x_new);
    change = norm(x_new - x_old);
    x_old = x_new;
end

% Initialize solution collections
u = [u0; v0];
u_collection = zeros(2*M,N);
x_collection = zeros(M,N);
H = zeros(N,1);

H(1) = Hamiltonian(u,x_new);
x_collection(:,1) = x_old;
u_collection(:,1) = u;
if doplot
   figure 
end
for i = 1:N-1
    if moving 
        % Adjust grid according to de Boor's algorithm
        x_new = deBoor(x_old,u(1:M));
        
        % Interpolate onto new grid, linear or cubic
        if interptype == 1
            u_int = linint_noncyclic(u,x_old,x_new);
        else
            u_int = zeros(size(u));
            u_int(1:M) = pchip(x_old,u(1:M),x_new);
            u_int(M+1:end) = pchip(x_old,u(M+1:end),x_new);
        end
    else
        u_int = u;
    end
    
    % Save grid in collection
    x_collection(:,i+1) = x_new;
    x_old = x_new;
    
    % Hamiltonian difference after interpolation
    diffH = H(i) - Hamiltonian(u_int,x_new);
    
    % Take one time step
    u = u_int +0.01;%*rand(size(u_new));
    f = @(z) F(z,u_int,diffH,dt,x_new,projtype);
    options = optimset('Display','off','TolFun',1e-14,'Jacobian','on');
    [u,~,~,output] = fsolve(f,u,options);
    
    % Calculate Hamiltonian at new time step and save this step's solution
    H(i+1) = Hamiltonian(u,x_new);
    u_collection(:,i+1) = u;
    
    % If plotting, plot solution
    if ~mod(i,2) && doplot
        t = tmin + i*dt;
        subplot(2,1,1)
        plot(x_new,u_collection(1:M,i+1))
        hold off
        axis([x_old(1) x_old(end) -1 7])
        title(['t = ' int2str(t)])
        subplot(2,1,2)
        plot(x_new,u_collection(M+1:end,i),'r')
        hold off
        axis([x_old(1) x_old(end) -2 14])
        title(['t = ' int2str(t)])
        if have_u_analytic % If comparing to analytic solution
            subplot(2,1,1)
            hold on
            u_an = u_analytic(x_new,t);
            plot(x_new,u_an)
            hold off
        end
        if have_v_analytic % If comparing to analytic solution
            subplot(2,1,2)
            hold on
            v_an = v_analytic(x_new,t);
            plot(x_new,v_an)
            hold off
        end
        pause(0.01)
    end
end

% If plotting, plot Hamiltonian error and mesh trajectories
if doplot
    t = tmin:dt:tmax;
    figure
    plot(t,H-H(1))
    title('Hamiltonian error')
    xlabel('t')
    
    figure
    colors = hsv(M+1);
    for i = 1:3:M
        
        plot(x_collection(i,:),t,'Color',colors(i,:))
        hold on
    end
    title('Mesh trajectories')
    ylabel('t')
    xlabel('x')
end

end

