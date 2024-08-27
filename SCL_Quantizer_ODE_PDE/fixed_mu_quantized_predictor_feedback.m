%Reset
close all
clear all 
clc
global  A B D U % Mbar1 mu0 Omega t0 T M Delta t_f

t_y= 0; % initial time t_0
Ts=10;% time limit 

% Parameters for the ODE (A,B) must be controllable rank(ctrb(A, B))
A = [-1, 1; 
     0, 1];
B = [0;
    1];
K = -place(A,B,[-1,-2]);

D=1; % The delay constant 
% Space discretization Omega=[0,h] :
L=D; % the right bound of the domain  
Nx= 100; % the number of points of space discretization (number of mesh points)
x= linspace(0,L,Nx); % the space discretization
dx= x(2) - x(1); % the space discretization step

% Time discretization :
t_f= t_y+Ts; % final time for integration on intervale [t_i,t_i+Ts]
Nt= 100; % the number of points of time discretization (time)
howfar= t_f/Nt;
t= t_y; % the initialization time

coef= 0.65; % coef depends on lambda_max to satisfy CFL condition
timestep= coef*dx;% Time discretization's step

% Initialization of the solutions by zeros
Yevol1= zeros(Nt,Nx); % X_1(t,x)
Yevol2= zeros(Nt,Nx); % X_2(t,x)

tout_y= zeros(Nt,1); %initialize the time

Y= InitialConditions(x); %  Y(t0,x)
%Y= zeros(3,Nx);  
Yevol1(1,:)= Y(1,:); %  X_1(0), ..., X_1(0) Nx times
Yevol2(1,:)= Y(2,:); %  X_2(0), ..., X_2(0) Nx times
Yevol3(1,:)= Y(3,:); %  u(0,x)

mu_t=.1;% a small value of switching variable
%mu_t=100;a large value of switching variable
M=2;
Delta=M/100;
Mbar1= 2;
t0=.17;
U_values = zeros(Nt,1); % Initialize array to store U(t)
%mu_values = zeros(Nt,1); % Initialize array to store mu(t)

tic
% Defining the structure of the original system
sol= setup(1,@DefineTransportEquation,t,x,Y,'LxF',[],@TransportBoundaryConditions); 
Integral = zeros(size(x));
        for j=1:Nx
            Integral(1,j)= K*expm(D*A)*expm(-x(1,j)*A)*B*mu_t*quantizer(Y(3,j),mu_t,M,Delta); 
        end
     qX = [mu_t*quantizer(Y(1,1),mu_t, M, Delta);mu_t*quantizer(Y(2,1),mu_t, M, Delta)];   
     U= K*expm(D*A)*qX + trapz(x,Integral);
    
 for m= 2:Nt  
    sol= hpde(sol,howfar,timestep); % The solver
    t_y= sol.t; % The current time 
    Y= sol.u; % The solutions of the original system       
    if t_y<=t0
        U=0;
    else 
        for j=1:Nx
            Integral(1,j)= K*expm(D*A)*expm(-x(1,j)*A)*B*mu_t*quantizer(Y(3,j),mu_t,M,Delta); 
        end
     qX = [mu_t*quantizer(Y(1,1),mu_t, M, Delta);mu_t*quantizer(Y(2,1),mu_t, M, Delta)];   
     U= K*expm(D*A)*qX + trapz(x,Integral);
    end

    %fprintf('ty %.3f U%.2f\t',t_y,U)
    tout_y(m)= t_y; % Storage of the solutions to use later to plot the solutions    
    Yevol1(m,1)= Y(1,1); % X1   
    Yevol2(m,1)= Y(2,1); % X2   
    Yevol3(m,:)= Y(3,:); % u(t,x)   
    %fprintf(' u(x,t) %.2f \t', Y(3,1))
    U_values(m) = U;
    %fprintf(' U%.2f \t', U)
    % mu_values(m) = mu_t; % Store mu(t)
end  
toc

% Initialize storage for norms of X and u
norm_X = zeros(Nt, 1);
sup_norm_u = zeros(Nt, 1);
for m = 1:Nt
    % Compute norms of X and sup norm of u at each time step
    norm_X(m) = norm([Yevol1(m,1); Yevol2(m,1)]);
    sup_norm_u(m) = max(abs(Yevol3(m,:)));
end

% Plot state evolution of ODE over time
Nfig= 0;    
%%%%%%%%%Figure 1 : The solution X_1(t)
Nfig= Nfig+1; figure(Nfig)
% Plot X_1(t)
subplot(2,2,1);
plot(tout_y, Yevol1(:,1),'k','LineWidth', 2); % plot X_1(t)
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$X_1(t)$', 'Interpreter', 'latex', 'FontSize', 16);
% Plot X_2(t)
subplot(2,2,2);
plot(tout_y, Yevol2(:,1),'k','LineWidth', 2); % plot X_2(t)
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$X_2(t)$', 'Interpreter', 'latex', 'FontSize', 16);
% Plot coupled solution of transport equation and ODE
subplot(2,2,[3,4]);
pstep = 3; % Increase the spacing between spatial grid points in the plot
tstep = 3; % Increase the spacing between time steps in the plot
[Y, Z] = meshgrid(x(1:pstep:end), howfar*(0:tstep:Nt-1));
U_plot = Yevol3(1:pstep:end, 1:tstep:end);
mesh(Y, Z, U_plot, 'LineWidth', 1,'edgecolor', 'black');
view(83,10);
xlabel('$x$', 'interpreter', 'latex', 'FontSize', 16);
ylabel('$t$', 'interpreter', 'latex', 'FontSize', 16);
zlabel('$u(x,t)$', 'interpreter', 'latex', 'FontSize', 16);
%%%%%%%%% Figure 2
Nfig = Nfig + 1;
figure(Nfig);
plot(tout_y, sup_norm_u+norm_X,'k','LineWidth', 3);
xlabel('$t$','Interpreter','latex','FontSize',16);
ylabel('$|X(t)|+\|u(t)\|_{\infty}$','Interpreter','latex','FontSize',16);


% Definition of the inital condition of ODE-PDE system
function Y = InitialConditions(x)
      Y(1,:)= 0*x + 10; %X_1(0)
      Y(2,:)= 0*x + 0; %X_2(0)
    Y(3,:)= 10; % u(0,x)%
end   
% Definition of the linear ODE-PDE system
function Y_t= DefineTransportEquation(t,x,Y,Y_x) 
    global A B
    Y_t(1,:)=  A(1,1)*Y(1,:) + A(1,2)*Y(2,:) + B(1,1)*Y(3,1); % ODE X_1(t)
    Y_t(2,:)=  A(2,1)*Y(1,:) + A(2,2)*Y(2,:) + B(2,1)*Y(3,1); % ODE X_2(t)

    Y_t(3,:)=  Y_x(3,:)+0*Y(3,:); % PDE u(t,x)
end
% end fucntion   

% Definition of the boundary condition of the ODE-PDE system
function [YL,YR]= TransportBoundaryConditions(t,YLex,YRex) 
    global U
    YL(1)= YLex(1); % No left boundary X_1(t)
    YR(1)= YRex(1); % No left boundary X_1(t)

    YL(2)= YLex(2); % No left boundary X_2(t)
    YR(2)= YRex(2); % No left boundary X_2(t)

    YL(3)= YLex(3); % No left boundary u(t,x)
    YR(3)= U;  % Control
end

% quantizer
function quantized_u = quantizer(u,mu_t, M, Delta)
    if u/mu_t >= M
        quantized_u= M;
    elseif u/mu_t <= -M
        quantized_u= -M;
    else
        quantized_u= Delta* floor(u/(Delta*mu_t) + 0.5);
    end
end