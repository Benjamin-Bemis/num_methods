%{

@author: Benjamin Bemis Ph.D Student, 
Advisor: Dr Juliano


Description:
AME 60614: Numerical Methods
Homework: 3
Due: 10/10/2024


%}

%%
clear all 
clc 
close all


%% Preperation of the workspace
fontsize = 16;


% set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultAxesFontSize',fontsize)
set(0,'DefaultLegendFontSize',fontsize)
colors  = ["#000000","#1b9e77","#d95f02","#7570b3","#0099FF"];


%% Problem 1




%% Problem 2

func =@(X) (100/(sqrt(X+0.01)))+1/((X-0.3)^2+0.001) - pi;
x = linspace(0,1,1000);

for i = 1:length(x)
    f(i) = (100/(sqrt(x(i)+0.01)))+1/((x(i)-0.3)^2+0.001) - pi;
end

figure
plot(x,f)
xlabel("$x$")
ylabel("$f(x)$")

%{
Functions to generate  



%}

n=[10:10:1e3];
syms X
F=int((100/(sqrt(X+0.01)))+1/((X-0.3)^2+0.001) - pi,0,1);
t2=double(F)*ones(size(n));
for i=1:length(n);
    x=linspace(0,1,n(i));
    f = (100./(sqrt(x+0.01)))+1./((x-0.3).^2+0.001) - pi;
    t(i)=trapz(x,f);
    simp(i) = simpson(func,0,1,n(i));
    t_end(i) = trapz_end(func,0,1,n(i));
    error(i)=abs(t(i)-t2(i));
    error_simp(i)=abs(simp(i)-t2(i));
    error_end(i)=abs(t_end(i)-t2(i));
end
figure
semilogy(n,error)
hold on 
semilogy(n,error_simp)
hold on 
semilogy(n,error_end)
xlabel("$n$")
ylabel("$I_{exact}-I_{num}$")
legend("Trap","Simpson","Trap Corrected")

%% Problem 3

func3 = @(X) (X+0.5).^(-2);
dfunc3 = diff(func3,X);
dfunc3_ex = [double(subs(dfunc3,1)),double(subs(dfunc3,5))];

ho = 0.5;
tolerance = 1e-5;

x = linspace(0,6,1e3);
func3_plot = func3(x);

figure
plot(x,func3_plot)


[value(1),h(1),error(1)] = rich_extrap(func3,1,ho,dfunc3_ex(1),tolerance);
[value(2),h(2),error(2)] = rich_extrap(func3,5,ho,dfunc3_ex(2),tolerance);

%% Functions

function [value] = simpson(f,a,b,n)

h=(b-a)/n;
if rem(n,2)==1
   fprintf('\n Given n was odd. Creating even sub-interval:'); 
   n=n+1;
end
so=0;
se=0;
for k=1:1:n-1
    x(k)=a+k*h;
    y(k)=f(x(k));
    if rem(k,2)==1
       so=so+y(k);%sum of odd terms
     else
       se=se+y(k); %sum of even terms
    end
end
% Formula:  (h/3)*[(y0+yn)+4*(y1+y3+y5+..odd term)+2*(y2+y4+y6+...even terms)]
value = double(h/3*(f(a)+f(b)+4*so+2*se));

end


function [value] = trapz_end(f,a,b,n)

h=(b-a)/n;

inner_sum = 0;
for k=1:1:n-1
    x(k)=a+k*h;
    y(k)=f(x(k));
    inner_sum = inner_sum+y(k);
end
syms X
df = diff(f,X);
df1 = double(subs(df,a));
df2 = double(subs(df,b));
value = double( h/2*(f(a)+f(b)+2*inner_sum) - h^2/12 *(df2-df1) );

end

function [df_rich,h,error] = rich_extrap(f,x,ho,exact,tolerance)
    h = ho;
    D = zeros(1, 4); %Store D(h), D(h/2), D(h/4), D(h/8)
    
    % Initial central difference approximations at different step sizes
    D(1) = (f(x + h) - f(x - h)) / (2 * h);        % D(h)
    D(2) = (f(x + h/2) - f(x - h/2)) / (h);        % D(h/2)
    D(3) = (f(x + h/4) - f(x - h/4)) / (h/2);      % D(h/4)
    D(4) = (f(x + h/8) - f(x - h/8)) / (h/4);      % D(h/8)
    
    % Apply Richardson extrapolation
    for i = 1:3
        for j = 1:(4-i)
            D(j) = (2^(i+1) * D(j+1) - D(j)) / (2^(i+1) - 1);
        end
    end
        df_rich = D(1);
    
    error = abs(df_rich - exact);
    h_current = h;
    
    while error > tolerance
        h_current = h_current / 2;
        
        % Recompute central differences
        D(1) = (f(x + h_current) - f(x - h_current)) / (2 * h_current);
        D(2) = (f(x + h_current/2) - f(x - h_current/2)) / h_current;
        D(3) = (f(x + h_current/4) - f(x - h_current/4)) / (h_current/2);
        D(4) = (f(x + h_current/8) - f(x - h_current/8)) / (h_current/4);
        
        % Apply Richardson extrapolation again for higher order
        for i = 1:3
            for j = 1:(4-i)
                D(j) = (2^(i+1) * D(j+1) - D(j)) / (2^(i+1) - 1);
            end
        end
        
        % Check the error between successive iterations
        error = abs(df_rich - exact);
        df_rich = D(1);  % Update the current best approximation
    end
h = h_current
end


function df_richardson = richardson_extrapolation_8th_order(f, x, h, tol)
    % Preallocate array to store Richardson approximations at different step sizes
    D = zeros(1, 4); % We will store D(h), D(h/2), D(h/4), D(h/8)
    
    % Initial central difference approximations at different step sizes
    D(1) = (f(x + h) - f(x - h)) / (2 * h);        % D(h)
    D(2) = (f(x + h/2) - f(x - h/2)) / (h);        % D(h/2)
    D(3) = (f(x + h/4) - f(x - h/4)) / (h/2);      % D(h/4)
    D(4) = (f(x + h/8) - f(x - h/8)) / (h/4);      % D(h/8)
    
    % Apply Richardson extrapolation to increase accuracy
    for i = 1:3
        for j = 1:(4-i)
            D(j) = (2^(i+1) * D(j+1) - D(j)) / (2^(i+1) - 1);
        end
    end
    
    % Final approximation is the most accurate Richardson extrapolated result
    df_richardson = D(1);
    
    % Optional: iterative refinement until the error is within the tolerance
    error = inf;
    h_current = h;
    
    while error > tol
        h_current = h_current / 2; % Halve the step size
        
        % Recompute central differences
        D(1) = (f(x + h_current) - f(x - h_current)) / (2 * h_current);
        D(2) = (f(x + h_current/2) - f(x - h_current/2)) / h_current;
        D(3) = (f(x + h_current/4) - f(x - h_current/4)) / (h_current/2);
        D(4) = (f(x + h_current/8) - f(x - h_current/8)) / (h_current/4);
        
        % Apply Richardson extrapolation again for higher order
        for i = 1:3
            for j = 1:(4-i)
                D(j) = (2^(i+1) * D(j+1) - D(j)) / (2^(i+1) - 1);
            end
        end
        
        % Check the error between successive iterations
        error = abs(df_richardson - D(1));
        df_richardson = D(1);  % Update the current best approximation
    end
end

