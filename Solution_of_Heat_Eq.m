%% This code is written to solve a parabolic PDE, heat equation numerically.
% As numeric method simple explicit method is used by following
% Chapra,S.C.,and Canale,R.P., Numerical Methods for Engineers, 7th ed.,...
% McGraw-Hill Education, 2015.

%% CASE 5
clear; clc; close all;

%Givens
L = 1;    %Length (m)
N = [19,9,19,19,19];   %Number of space grids

alpha = 0.2;    % Coefficient of in	the	governing equation
time = 1.5;        % Time (s)
CFL = [0.25;0.5;0.5;0;2];

% Analytical Solution
syms x_2 t_2

u = @(x_2,t_2) 100*exp(-alpha*pi^2*t_2)*sin(pi*x_2);

for j = 1:length(N)
    x = 0;
    t = 0;
    dx = L/(N(j)-1);   % Equal space difference between nodes (m)
    x = 0:dx:L;
    
    if j == 4   % For case 4
        dt = alpha * dx;    % Largest possible time step for stability (s)
        CFL(j) = alpha * dt / dx^2; % Finding CFL with formulation
        t = 0:dt:time;
    else
        
        dt = (CFL(j)*dx^2)/alpha;  % Time difference(s)
        t = 0:dt:time;
    end
    
    M(j) = length(t);     %Number of time grids
    
    T_e = zeros(N(j),M(j));
    
    %Boundary Condition
    T_e(1,:) = 0;
    T_e(N(j),:) = 0;
    
    % %Initial Condition
    for i = 2:N(j)-1
        T_e(i,1) = 100 * sin(pi*x(i)/L);
    end
    
    T_i = T_e;
    T_i = implicit(CFL(j),N(j),M(j),T_i);
    
    %Main Formula
    for k = 1:M(j)-1
        for i = 2:N(j)-1
            T_e(i,k+1) = T_e(i,k) + CFL(j) * (T_e(i+1,k) - 2* T_e(i,k) + T_e(i-1,k));
        end
    end
    
    Case(j).x = x';
    Case(j).t = t';
    Case(j).T_e = T_e';
    Case(j).T_i = T_i';
    
    T_analytical = zeros(N(j),M(j));
    
    for i = 1:N(j)
        for k = 1:M(j)
            T_analytical(i,k) = u(Case(j).x(i),Case(j).t(k));
        end
    end
    
    Case(j).T_an = T_analytical';
    
    error = abs((T_e-T_analytical) ./ (T_e)) *100;  % Percentage absulate error
    error(isinf(error)) = 0;
    error(isnan(error)) = 0;
    
    Case(j).error_e = error';
    
    error = abs((T_i-T_analytical) ./ (T_i)) *100;  % Percentage absulate error
    error(isinf(error)) = 0;
    error(isnan(error)) = 0;
    
    Case(j).error_i = error';
end

%RMSE Calculation
for j = 1:length(N)
    RMSE(j).a = sqrt(sum(abs(Case(j).T_an - Case(j).T_e).^2 )) / M(j);
end

figure(1)
fplot(u(x_2,0),[0,1],'b')
hold on

plot(Case(1).x,Case(1).T_e(1,:),'--.','MarkerSize',10)
hold on

plot(Case(2).x,Case(2).T_e(1,:),'--.','MarkerSize',10)
hold on

plot(Case(3).x,Case(3).T_e(1,:),'--.','MarkerSize',10)
hold on

% plot(Case(4).x,Case(4).T_e(1,:),'c--^','MarkerSize',10)
% hold on

fplot(u(x_2,1.5),[0,1],'b')
hold on

plot(Case(1).x,Case(1).T_e(M(1),:),'--.','MarkerSize',10)
hold on

plot(Case(2).x,Case(2).T_e(M(2),:),'--.','MarkerSize',10)
hold on

plot(Case(3).x,Case(3).T_e(M(3),:),'--.','MarkerSize',10)
hold on

% plot(Case(4).x,Case(4).T_e(M(4),:),'c--^','MarkerSize',10)
% hold on
%
grid on
legend('Analytical t=0',...
    'Case1 t=0',...
    'Case2 t=0',...
    'Case3 t=0',...
    'Analytical t =1.5',...
    'Case1 t=1.5',...
    'Case2 t=1.5',...
    'Case3 t=1.5',...
    'Interpreter','Latex')
xlabel('x $(m)$','Interpreter','Latex')
ylabel('Temperature $(^{\circ}C)$','Interpreter','Latex')
title('Solution with Simple Explicit Method','Interpreter','Latex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
figure(2)
fplot(u(x_2,0),[0,1],'b')
hold on

plot(Case(1).x,Case(1).T_i(1,:),'--.','MarkerSize',10)
hold on

plot(Case(2).x,Case(2).T_i(1,:),'--.','MarkerSize',10)
hold on

plot(Case(3).x,Case(3).T_i(1,:),'--.','MarkerSize',10)
hold on

fplot(u(x_2,1.5),[0,1],'b')
hold on

plot(Case(1).x,Case(1).T_i(M(1),:),'--.','MarkerSize',10)
hold on

plot(Case(2).x,Case(2).T_i(M(2),:),'--.','MarkerSize',10)
hold on

plot(Case(3).x,Case(3).T_i(M(3),:),'--.','MarkerSize',10)
hold on
grid on
legend('Analytical t=0',...
    'Case1 t=0',...
    'Case2 t=0',...
    'Case3 t=0',...
    'Analytical t =1.5',...
    'Case1 t=1.5',...
    'Case2 t=1.5',...
    'Case3 t=1.5',...
    'Interpreter','Latex')
xlabel('x $(m)$','Interpreter','Latex')
ylabel('Temperature $(^{\circ}C)$','Interpreter','Latex')
title('Solution with Simple Implicit Method','Interpreter','Latex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
subplot(2,1,1)
plot(Case(4).x,Case(4).T_e(M(4),:),'--o','MarkerSize',10)
hold on
plot(Case(4).x,Case(4).T_e(1,:),'--o','MarkerSize',10)
hold on
plot(Case(5).x,Case(5).T_e(M(5),:),'--o','MarkerSize',10)
hold on
plot(Case(5).x,Case(5).T_e(1,:),'--o','MarkerSize',10)
hold on
fplot(u(x_2,1.5),[0,1],'b')
hold on
fplot(u(x_2,0),[0,1],'b')
grid on

legend('Case4 t=1.5','Case4 t=0',...
    'Case5 t=1.5','Case5 t=0', 'Analytical t=1.5', 'Analytical t=0' )
xlabel('x $(m)$','Interpreter','Latex')
ylabel('Temperature $(^{\circ}C)$','Interpreter','Latex')
title('Stability Comparison of Methods','Interpreter','Latex')

subplot(2,1,2)
plot(Case(4).x,Case(4).T_i(M(4),:),'--o','MarkerSize',10)
hold on
plot(Case(4).x,Case(4).T_i(1,:),'--o','MarkerSize',10)
hold on
plot(Case(5).x,Case(5).T_i(M(5),:),'--o','MarkerSize',10)
hold on
plot(Case(5).x,Case(5).T_i(1,:),'--o','MarkerSize',10)
hold on
fplot(u(x_2,1.5),[0,1],'b')
hold on
fplot(u(x_2,0),[0,1],'b')
grid on

legend('Case4 t=1.5','Case4 t=0',...
    'Case5 t=1.5','Case5 t=0', 'Analytical t=1.5', 'Analytical t=0' )
xlabel('x $(m)$','Interpreter','Latex')
ylabel('Temperature $(^{\circ}C)$','Interpreter','Latex')
