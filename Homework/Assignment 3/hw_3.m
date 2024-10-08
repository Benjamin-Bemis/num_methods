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
    error(i)=abs(t2(i)-t(i));
    error_simp(i)=abs(t2(i)-simp(i));
    error_end(i)=abs(t2(i)-t_end(i));
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