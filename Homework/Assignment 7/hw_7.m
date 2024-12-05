%{

@author: Benjamin Bemis Ph.D Student, 
Advisor: Dr Juliano


Description:
AME 60614: Numerical Methods
Homework: 7
Due: 12/5/2024


%}

%% Preperation of the workspace
clear all 
clc 
close all
fontsize = 16;


% set(0,'DefaultFigureWindowStyle','default')
set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultAxesFontSize',fontsize)
set(0,'DefaultLegendFontSize',fontsize)
colors  = ["#000000","#1b9e77","#d95f02","#7570b3","#0099FF"]';

%% Setting data paths
% Make sure to update this for the machine that you are working on. (Maybe, This should now run on any machine without change. 7/24/24)
% Change the current folder to the folder of this m-file.
if(~isdeployed)
  cd(fileparts(matlab.desktop.editor.getActiveFilename));
end

addpath(cd)
% cd ..; % Moving up a directory (from processing_code)
basepath = cd; % Pulling the current directory


imagepath = [basepath filesep 'images' filesep]; 
mkdir(imagepath);

   


%% Problem 9 Chapter 5

L = 1;                   % Length of the domain (0 ≤ x ≤ 1)
Nx = 51;                 % Number of spatial points
dx = L / (Nx - 1);       % Spatial step size
x = linspace(0, L, Nx);  % Spatial grid
u = 0.08;                % Convection velocity
dt_values = [0.0625, 0.03, 0.001];  % Try different dt values
t_end = 8;               % End time
Nt = round(t_end / dt_values(1));   % Number of time steps for largest dt

% Exact solution function
exact_solution = @(x, t, u) (x - u * t >= 0 & x - u * t <= 0.2) .* (1 - (10 * (x - u * t) - 1).^2);

% Initial condition
T = zeros(Nx, 1);
for i = 1:Nx
    if x(i) <= 0.2
        T(i) = 1 - (10 * x(i) - 1)^2;
    end
end

% Plot initial condition
figure;
plot(x, T, '-', 'LineWidth', 1.5, 'Color', colors(1));
hold on;
title('Pure Convection: Initial and Exact Solutions');
xlabel('x'); ylabel('T(x,t)');
legend('Initial Condition');

% Loop through different time steps for comparison
for k = 1:length(dt_values)
    dt = dt_values(k);
    T_num = T;            % Numerical solution for current dt
    Nt = round(t_end / dt); % Number of time steps

    % Stability check for explicit Euler (should satisfy CFL condition u*dt/dx <= 1)
    if u * dt / dx > 1
        warning('CFL condition violated: u*dt/dx > 1. Results may be unstable.');
    end

    % Time marching using Explicit Euler scheme
    for n = 1:Nt
        T_new = T_num;  % Temporary array to store updated values

        % Update T using the central difference for spatial derivative
        for i = 2:Nx-1
            T_new(i) = T_num(i) - (u * dt / (2 * dx)) * (T_num(i+1) - T_num(i-1));
        end

        % Update boundary conditions (Dirichlet T(0, t) = 0, T(1, t) = 0)
        T_new(1) = 0;
        T_new(Nx) = 0;

        % Update the solution for the next time step
        T_num = T_new;
    end

    for t_plot = [0, 4, 8]
        T_exact = exact_solution(x, t_plot, u);  % Exact solution
        figure
        plot(x, T_exact, '-', 'LineWidth', 1.5, 'Color', colors(1))
        hold on
        plot(x, T_num, '--', 'LineWidth', 1.5, 'Color',colors(2))        
        title(['Solution Comparison at t = ', num2str(t_plot), ' for dt = ', num2str(dt)]);
        xlabel('x'); ylabel('T(x,t)');
        legend('Exact Solution','Numerical Solution');
    end
end



%% Functions

function [root,err] = Nraph(f,initGuess,tol)
%Nraph solves for the root nearest the initial guess using the Newt-Raphson
%method.

% INPUTS:
% f: is a function handle. 
% initGuess: is the inital guess for the root.
% tol: desired tolerance.

% OUTPUTS:
% root: nearest root to initial guess.
% err: error in the solution of the root.

x = initGuess;
x(2) = initGuess+1;
counter = 2;
while abs(f(x(end))) >= tol

x(counter+1) = x(counter) - f(x(counter)) * (x(counter)-x(counter-1))/(f(x(counter))-f(x(counter-1)));

counter = counter+1;
end
root = x(counter);
err =  abs(f(root));

end

function [x, y1, y2] = shoot(f,xRange, BC, h, dx_guess ,tol)

% x = xRange(1):h:xRange(2);

dx = dx_guess;
dx(2) = dx_guess+.5;

% 
y1 = inf;
% [x, y1, y2] = RK4_2(f, BC(1), dx_guess, xRange(1), xRange(2), h);

counter = 1;
while abs(y1(end)- BC(2)) >= tol

    if counter == 1
        [x, y1, y2] = RK4_2(f, BC(1), dx_guess, xRange(1), xRange(2), h);
        disp(["counter = ", string(counter)])

    elseif counter == 2
        [x, y1, y2] = RK4_2(f, BC(1), dx_guess+.5, xRange(1), xRange(2), h);
        disp(["counter = ", string(counter)])
         dx(counter+1) = dx(counter) - (y1(end)-BC(2)) * (dx(counter) - dx(counter-1))/ (y1(end)-y_prev(counter));
    elseif counter > 2
        [x, y1, y2] = RK4_2(f, BC(1), dx(counter), xRange(1), xRange(2), h);
        dx(counter+1) = dx(counter) - (y1(end)-BC(2)) * (dx(counter) - dx(counter-1))/ (y1(end)-y_prev(counter));

    end

y_prev(counter+1) = y1(end);
err = abs(y1(end)- BC(2))
counter = counter+1;
end
dx(end-1)
end

function [x, y1, y2] = shootdf(f,xRange, BC, h, dx_guess ,tol)

% x = xRange(1):h:xRange(2);

dx = dx_guess;
dx(2) = dx_guess+.5;

% 
y2 = inf;
% [x, y1, y2] = RK4_2(f, BC(1), dx_guess, xRange(1), xRange(2), h);

counter = 1;
while abs(y2(end)- BC(2)) >= tol

    if counter == 1
        [x, y1, y2] = RK4_2(f, BC(1), dx_guess, xRange(1), xRange(2), h);
        disp(["counter = ", string(counter)])

    elseif counter == 2
        [x, y1, y2] = RK4_2(f, BC(1), dx_guess+.5, xRange(1), xRange(2), h);
        disp(["counter = ", string(counter)])
         dx(counter+1) = dx(counter) - (y2(end)-BC(2)) * (dx(counter) - dx(counter-1))/ (y2(end)-y_prev(counter));
    elseif counter > 2
        [x, y1, y2] = RK4_2(f, BC(1), dx(counter), xRange(1), xRange(2), h);
        dx(counter+1) = dx(counter) - (y2(end)-BC(2)) * (dx(counter) - dx(counter-1))/ (y2(end)-y_prev(counter));

    end

y_prev(counter+1) = y2(end);
err = abs(y2(end)- BC(2))
counter = counter+1;
end
dx(end-1)
end

function [t, y1, y2] = RK4_2(f, y0, v0, t0, tf, h)

    % RK4_2 works for single equation odes
 
    % Inputs:
    %   f  - Function handle for y'' = f(t, y, y')
    %   y0 - Initial condition for y (y(t0) = y0)
    %   v0 - Initial condition for y' (y'(t0) = v0)
    %   t0 - Initial time
    %   tf - Final time
    %   h  - Step size
    %
    % Outputs:
    %   t  - Array of time steps
    %   y1 - Array of solution values for y at each time step
    %   y2 - Array of solution values for y' at each time step

    % Define the time vector from t0 to tf with step size h
    t = t0:h:tf;
    N = length(t); % Number of time steps
    y1 = zeros(1, N); % Preallocate y1 for y
    y2 = zeros(1, N); % Preallocate y2 for y'

    % Set the initial conditions
    y1(1) = y0;
    y2(1) = v0;

    % Apply the 4th-order Runge-Kutta method
    for n = 1:N-1
        % Calculate k1 values
        k1_y1 = y2(n);
        k1_y2 = f(t(n), y1(n), y2(n));

        % Calculate k2 values
        k2_y1 = y2(n) + h/2 * k1_y2;
        k2_y2 = f(t(n) + h/2, y1(n) + h/2 * k1_y1, y2(n) + h/2 * k1_y2);

        % Calculate k3 values
        k3_y1 = y2(n) + h/2 * k2_y2;
        k3_y2 = f(t(n) + h/2, y1(n) + h/2 * k2_y1, y2(n) + h/2 * k2_y2);

        % Calculate k4 values
        k4_y1 = y2(n) + h * k3_y2;
        k4_y2 = f(t(n) + h, y1(n) + h * k3_y1, y2(n) + h * k3_y2);

        % Update y1 and y2 using weighted average of slopes
        y1(n+1) = y1(n) + (h/6) * (k1_y1 + 2*k2_y1 + 2*k3_y1 + k4_y1);
        y2(n+1) = y2(n) + (h/6) * (k1_y2 + 2*k2_y2 + 2*k3_y2 + k4_y2);
    end
end

function [t, Y] = RK4(f, Y0, t0, tf, h)
    % RK4 - 4th-order Runge-Kutta method for systems of equations.
    %
    % Inputs:
    %   f  - Function handle for the system of equations, f(t, Y)
    %        Y is a column vector, and f should return a column vector.
    %   Y0 - Initial conditions as a column vector (Nx1, where N is the number of equations)
    %   t0 - Initial time
    %   tf - Final time
    %   h  - Step size
    %
    % Outputs:
    %   t  - Array of time steps
    %   Y  - Solution matrix (NxM, where N is the number of equations, M is the number of time steps)

    % Define the time vector from t0 to tf with step size h
    t = t0:h:tf;
    N = length(t);       % Number of time steps
    num_eqns = length(Y0);  % Number of equations in the system

    % Preallocate the solution matrix Y
    Y = zeros(num_eqns, N);
    
    % Set the initial conditions
    Y(:, 1) = Y0;

    % Apply the 4th-order Runge-Kutta method for each time step
    for n = 1:N-1
        % Calculate k1 values
        k1 = f(t(n), Y(:, n));

        % Calculate k2 values
        k2 = f(t(n) + h/2, Y(:, n) + h/2 * k1);

        % Calculate k3 values
        k3 = f(t(n) + h/2, Y(:, n) + h/2 * k2);

        % Calculate k4 values
        k4 = f(t(n) + h, Y(:, n) + h * k3);

        % Update Y using the weighted average of slopes
        Y(:, n+1) = Y(:, n) + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
    end
end

function [t, Y] = trapezoidal_s(f, Y0, t0, tf, h)
    % trapezoidal_s - Linearized Trapezoidal method for systems of ODEs.
    %
    % Inputs:
    %   f  - Function handle for the system of equations, f(t, Y)
    %        Y is a column vector, and f should return a column vector.
    %   Y0 - Initial conditions as a column vector (Nx1, where N is the number of equations)
    %   t0 - Initial time
    %   tf - Final time
    %   h  - Step size
    %
    % Outputs:
    %   t  - Array of time steps
    %   Y  - Solution matrix (NxM, where N is the number of equations, M is the number of time steps)

    % Define the time vector from t0 to tf with step size h
    t = t0:h:tf;
    N = length(t);       % Number of time steps
    num_eqns = length(Y0);  % Number of equations in the system

    % Preallocate the solution matrix Y
    Y = zeros(num_eqns, N);
    
    % Set the initial conditions
    Y(:, 1) = Y0;

    % Apply the Linearized Trapezoidal method for each time step
    for n = 1:N-1
        % Predictor step (Euler's method)
        Y_star = Y(:, n) + h * f(t(n), Y(:, n));
        
        % Corrector step
        Y(:, n+1) = Y(:, n) + (h/2) * (f(t(n), Y(:, n)) + f(t(n+1), Y_star));
    end
end
