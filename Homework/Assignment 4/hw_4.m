%{

@author: Benjamin Bemis Ph.D Student, 
Advisor: Dr Juliano


Description:
AME 60614: Numerical Methods
Homework: 4
Due: 10/31/2024


%}

%% Preperation of the workspace
clear all 
clc 
close all
fontsize = 16;


% set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultAxesFontSize',fontsize)
set(0,'DefaultLegendFontSize',fontsize)
colors  = ["#000000","#1b9e77","#d95f02","#7570b3","#0099FF"]';

%% Problem 2 in Chapter 4

f = @(t,y) (-0.2*y) - (2*cos(2*t)*y^2);
y0 = 1;
t0 = 0;
tf = 7;
h = [0.2 0.05 0.025 0.006];

figure
for n = 1:size(h,2)
       [t,y] = explicitEuler(f, y0, t0, tf, h(n));
       plot(t,y, "LineWidth",2, 'DisplayName',strcat("h = ",string(h(n))) , color=colors(n,:))
       hold on
end
xlabel('$t$ (sec)')
ylabel('$v$')
grid on
xlim([t0 tf])
ylim([0 1.4])
legend(Location="best",Interpreter="latex")


%% Problem 6 in Chapter 4

f = @(t,y) -((3*t)/(1+t))*y - (2*(1+t)^3*exp(-t));
y0 = 1;
t0 = 0;
tf = 15;
h = [0.2 0.8 1.1];





%% Problem 8 in Chapter 4





%% Functions

function [t, y] = explicitEuler(f, y0, t0, tf, h)
    % explicitEuler solves an ODE using the explicit Euler method.
    %
    % Inputs:
    %   f  - Function handle for dy/dt = f(t, y)
    %   y0 - Initial condition (value of y at t = t0)
    %   t0 - Initial time
    %   tf - Final time
    %   h  - Step size
    %
    % Outputs:
    %   t - Array of time steps
    %   y - Array of solution values at each time step

    % Define the time vector from t0 to tf with step size h
    t = t0:h:tf;
    N = length(t); % Number of time steps
    y = zeros(1, N); % Preallocate y for speed

    % Set the initial condition
    y(1) = y0;

    % Apply the explicit Euler method
    for n = 1:N-1
        y(n+1) = y(n) + h * f(t(n), y(n));
    end
end

function [t, y] = implicitEuler(f, y0, t0, tf, h)
    % implicitEuler solves an ODE using the implicit Euler method.
    %
    % Inputs:
    %   f  - Function handle for dy/dt = f(t, y)
    %   y0 - Initial condition (value of y at t = t0)
    %   t0 - Initial time
    %   tf - Final time
    %   h  - Step size
    %
    % Outputs:
    %   t - Array of time steps
    %   y - Array of solution values at each time step

    % Define the time vector from t0 to tf with step size h
    t = t0:h:tf;
    N = length(t); % Number of time steps
    y = zeros(1, N); % Preallocate y for speed

    % Set the initial condition
    y(1) = y0;

    % Options for fsolve to increase accuracy and ensure convergence
    options = optimoptions('fsolve', 'Display', 'off');

    % Apply the implicit Euler method
    for n = 1:N-1
        % Define the function for the nonlinear equation at each step
        g = @(ynext) ynext - y(n) - h * f(t(n+1), ynext);
        
        % Use fsolve to solve for y(n+1)
        y(n+1) = fsolve(g, y(n), options);
    end
end

function [t, y] = trapMethod(f, y0, t0, tf, h)
    % trapMethod solves an ODE using the trapezoidal method.
    %
    % Inputs:
    %   f  - Function handle for dy/dt = f(t, y)
    %   y0 - Initial condition (value of y at t = t0)
    %   t0 - Initial time
    %   tf - Final time
    %   h  - Step size
    %
    % Outputs:
    %   t - Array of time steps
    %   y - Array of solution values at each time step

    % Define the time vector from t0 to tf with step size h
    t = t0:h:tf;
    N = length(t); % Number of time steps
    y = zeros(1, N); % Preallocate y for speed

    % Set the initial condition
    y(1) = y0;

    % Options for fsolve to increase accuracy and ensure convergence
    options = optimoptions('fsolve', 'Display', 'off');

    % Apply the trapezoidal method
    for n = 1:N-1
        % Define the function for the nonlinear equation at each step
        g = @(ynext) ynext - y(n) - (h/2) * (f(t(n), y(n)) + f(t(n+1), ynext));
        
        % Use fsolve to solve for y(n+1)
        y(n+1) = fsolve(g, y(n), options);
    end
end

function [t, y] = RK2(f, y0, t0, tf, h)
    % RK2 solves an ODE using the 2nd-order Runge-Kutta method.
    %
    % Inputs:
    %   f  - Function handle for dy/dt = f(t, y)
    %   y0 - Initial condition (value of y at t = t0)
    %   t0 - Initial time
    %   tf - Final time
    %   h  - Step size
    %
    % Outputs:
    %   t - Array of time steps
    %   y - Array of solution values at each time step

    % Define the time vector from t0 to tf with step size h
    t = t0:h:tf;
    N = length(t); % Number of time steps
    y = zeros(1, N); % Preallocate y for speed

    % Set the initial condition
    y(1) = y0;

    % Apply the 2nd-order Runge-Kutta method
    for n = 1:N-1
        k1 = f(t(n), y(n));
        k2 = f(t(n) + h/2, y(n) + h/2 * k1);
        y(n+1) = y(n) + h * k2;
    end
end

function [t, y] = RK4(f, y0, t0, tf, h)
    % RK4 solves an ODE using the 4th-order Runge-Kutta method.
    %
    % Inputs:
    %   f  - Function handle for dy/dt = f(t, y)
    %   y0 - Initial condition (value of y at t = t0)
    %   t0 - Initial time
    %   tf - Final time
    %   h  - Step size
    %
    % Outputs:
    %   t - Array of time steps
    %   y - Array of solution values at each time step

    % Define the time vector from t0 to tf with step size h
    t = t0:h:tf;
    N = length(t); % Number of time steps
    y = zeros(1, N); % Preallocate y for speed

    % Set the initial condition
    y(1) = y0;

    % Apply the 4th-order Runge-Kutta method
    for n = 1:N-1
        k1 = f(t(n), y(n));
        k2 = f(t(n) + h/2, y(n) + h/2 * k1);
        k3 = f(t(n) + h/2, y(n) + h/2 * k2);
        k4 = f(t(n) + h, y(n) + h * k3);
        y(n+1) = y(n) + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
    end
end
