%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Simultaneous compensation of input delay and state quantization for  %
%       linear systems via switched predictor feedback                    %
%                                                            April 2024   %
%                                                                         %
%                                                                         %
%             Florent Koudohode, Nikolaos Bekiaris-Liberis                %                                                     
%                                                                         %
%                                                                         %
%                                                                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reset
close all
clear all 
clc
global  A B D U  Mbar1 mu0 Omega t0 T M Delta t_f tau

t_y= 0; % initial time t_0
Ts= 15;% time limit 

% Parameters for the ODE (A,B) must be controllable rank(ctrb(A, B))
A = [-1, 1; 0, 1];
B = [0; 1];
K = -place(A,B,[-1,-2]);

D=1; % The delay constant 
L=D; % the right bound of the domain  
Nx= 100; % the number of points of space discretization (number of mesh points)
x= linspace(0,L,Nx); % the space discretization
dx= x(2) - x(1); % the space discretization step

% Time discretization :
t_f= t_y+Ts; % final time for integration on intervale [t_i,t_i+Ts]
Nt= 500; % the number of points of time discretization (time)
howfar= t_f/Nt;
t= t_y; % the initialization time

coef= 0.6; %coef depends on lambda_max to satisfy CFL condition
timestep= coef*dx; %Time discretization's step

% Initialization of the solutions by zeros
Yevol1= zeros(Nt,Nx); % X_1(t,x)
Yevol2= zeros(Nt,Nx); % X_2(t,x)

tout_y= zeros(Nt,1); %initialize the time

Y= InitialConditions(x); %  Y(t0,x)

Yevol1(1,:)= Y(1,:); %  X_1(0), ..., X_1(0) Nx times
Yevol2(1,:)= Y(2,:); %  X_2(0), ..., X_2(0) Nx times
Yevol3(1,:)= Y(3,:); %  u(0,x)

%%% Norm using the backstepping kernel
M1=norm(K,'inf')*exp(norm(A,'inf')*D)*max(1,D*norm(B,'inf'))+1;
M2=1/(norm(K,'inf')*exp(norm(A+B*K,'inf')*D)*max(1,D*norm(B,'inf'))+1);
% % Calculate M3
M3=norm(K,'inf')*exp(norm(A,'inf')*D)*(1+norm(B,'inf')*D);
% % Calculate barM1
Mbar1= 1+D*norm(B,'inf');
eigenvalue_ABK=eig(A+B*K);
% %Parameters to compute M0 Mbar Omega T
% sigma=-max(eigenvalue_ABK); % Msigma and sigma are such that (24) holds: |exp(A+B*K)t|<=Msigma exp(-sigma) 
% Msigma=norm(A+B*K); 
 lambda=8; %lambda for the small-gain condition  (25)
%lambda=12;violate the small-gain condition
%
% epsilon=.5; % for fading memory
% nu=sigma+1; % for ISS transport
% phi=((1+epsilon)/(1+lambda))*exp(D*(nu+1));
% varphi1=((1+epsilon)/(1-phi))*phi*norm(B)*(Msigma/sigma);
% C0= max(exp(D*nu+D)/(1+epsilon)/(1-varphi1),phi*Msigma/((1-phi)*(1-varphi1)));
% C1=max((Msigma/(1-varphi1)),((1+epsilon)*Msigma*exp(D*nu+D)*norm(B,'inf'))/((1-phi)*(1-varphi1)*sigma));
% M0=C0+C1;
Mbar=.6;% Mbar=M2/(M1*(1+M0));
% Choose Delta and M are such that condition (28) holds
 M=2;% Range of the quantizer
 Delta =M/100;%Quantization error

% % Calculate Omega <1

Omega=.63;% Omega=((1+lambda)*(1+M0)*Delta*M3)/(M2*M);
% delta=11; %delta is in (0,sigma)
% % Compute T
% T=-log(Omega/(1+M0))/delta;
T=2;
mu0 = 1; % Initial value of mu
tau = 1; % tau arbitrary 
% %%%%%%%%%%%%
% Calculate t0 
normX0 = norm([Yevol1(1,1); Yevol2(1,1)],'inf'); % Norm of X0
normu0= max(abs(Yevol3(1,1))); % Infinity norm of u0
t0=0;%t0 = (1 / norm(A,'inf')) * log((1 / mu0) * (normX0 + normu0) / (M * Mbar - 2 * Delta));


U_values = zeros(Nt,1); % Initialize array to store the control U(t)
mu_values = zeros(Nt,1); % Initialize array to store the switching parameter mu(t)
mu_values(1)=NaN;

% Defining the structure of the original system
sol= setup(1,@DefineTransportEquation,t,x,Y,'LxF',[],@TransportBoundaryConditions); 
tic
Integral = zeros(size(x));
if  t_y<=t0
            U = 0;
end % 
for m= 2:Nt  
    sol= hpde(sol,howfar,timestep); %The solver
    t_y= sol.t; %The current time 
    Y= sol.u; %The solutions of the original system in the current time
    mu_t=mu(t_y);
 if  t_y<=t0
            U = 0;
 else
    for j=1:Nx
        Integral(1,j)= K*expm(D*A)*expm(-x(1,j)*A)*B*mu_t*quantizer(Y(3,j),mu_t, M, Delta); 
    end
    quantizerX=[mu_t*quantizer(Y(1,1),mu_t, M, Delta);mu_t*quantizer(Y(2,1),mu_t, M, Delta)];    
    U= K*expm(D*A)*quantizerX+ trapz(x,Integral); % U(t_m)
 end
    tout_y(m)= t_y; % Storage of the solutions to use later to plot the solutions    
    Yevol1(m,1)= Y(1,1); % X1   
    Yevol2(m,1)= Y(2,1); % X2   
    Yevol3(m,:)= Y(3,:); % u(t,x)  
    U_values(m) = U;
    mu_values(m) = mu_t; % Store mu(t)
end  
toc

% Initialize storage for norms of X and u
norm_X = zeros(Nt, 1);
sup_norm_u = zeros(Nt, 1);
for m = 1:Nt
    % Compute norms of X and sup norm of u at each time step
    norm_X(m) = norm([Yevol1(m,1); Yevol2(m,1)],'inf');
    sup_norm_u(m) = max(abs(Yevol3(m,:)));
end


% Plot state evolution over time
Nfig= 0;    
%%%%%%%%%%%%%% Figure 1 : Solution ([X_1(t),X_2(t)],u(x,t))
Nfig= Nfig+1; figure(Nfig)
% Plot X_1(t)
subplot(2,2,1);
plot(tout_y, Yevol1(:,1),'k','LineWidth', 3 ); % plot X_1(t)
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$X_1(t)$', 'Interpreter', 'latex', 'FontSize', 20);
% Plot X_2(t)
subplot(2,2,2);
plot(tout_y, Yevol2(:,1),'k','LineWidth', 3); % plot X_2(t)
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$X_2(t)$', 'Interpreter', 'latex', 'FontSize', 20);
% Plot coupled solution of transport equation and ODE
subplot(2,2,[3,4]);
pstep = 3; % Increase the spacing between spatial grid points in the plot
tstep = 3; % Increase the spacing between time steps in the plot
[Y, Z] = meshgrid(x(1:pstep:end), howfar*(0:tstep:Nt-1));
U_plot = Yevol3(1:pstep:end, 1:tstep:end);
mesh(Y, Z, U_plot, 'LineWidth', 1,'edgecolor', 'black');
view(83,10);
xlabel('$x$', 'interpreter', 'latex', 'FontSize', 20);
ylabel('$t$', 'interpreter', 'latex', 'FontSize', 20);
zlabel('$u(x,t)$', 'interpreter', 'latex', 'FontSize', 20);
 %%%%%%%%%%%%% Figure 2
 Nfig = Nfig + 1;
 figure(Nfig);
 plot(tout_y, sup_norm_u+norm_X,'k',tout_y, M*Mbar*mu_values,'k--','LineWidth', 3);
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
%Switching parameter
function mu_t = mu(t)
global  A  Mbar1 mu0 Omega t0 T t_f tau; 
       if t<=t0 
          j=1;
          while (j-1)*tau > t && t>j*tau && j<=floor(t0/tau)+1
          j=j+1;
          end 
          mu_t=Mbar1*exp(2 * norm(A) * tau * j) * mu0;
      elseif t0<t && t<= t0+T           
          mu_t=mu(t0);
      else        
           i=2;
            test=t0+(i-1)*T<t && t<=t0+i*T;
          while test==0 && i<=floor((t_f-t0)/T)
          i=i+1;
          test=t0+(i-1)*T<t && t<=t0+i*T;
          end
          mu_t=Omega*mu(t0+(i-1)*T);
      end
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