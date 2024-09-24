%{

@author: Benjamin Bemis Ph.D Student, 
Advisor: Dr Juliano


Description:
AME 60614: Numerical Methods
Homework: 2
Due: 9/24/2024


%}

%% Preparation of the Workspace

clear all
clc
close all

%% Preperation of Figures 

fontsize = 16;
set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultAxesFontSize',fontsize)
set(0,'DefaultLegendFontSize',fontsize)
colors  = ["#000000","#1b9e77","#d95f02","#7570b3","#0099FF","#FF0000"];


%% problem 1


%% problem 2
syms h
A = [ 1 1 1 1; -h 0 h 2*h; h^2/2 0 h^2/2 2*h^2; -h^3/6 0 h^3/6 (2*h)^3/6];
B = [0;-1;0;0];

X = linsolve(A,B)

%% problem 3

% part a
syms h a
A2 = [1 1 1 1; 0 h 2*h 3*h; 0 h^2/2 2*h^2 (3*h)^2/2;0 h^3/6 (2*h)^3/6 (3*h)^3/6];
B2 = [0;-1-a;-a*h;-a*h^2/2];

X2 = linsolve(A2,B2)


% part b

syms h
A3 = [0 1 1 1 1; 1 0 h 2*h 3*h; h 0 h^2/2 2*h^2 (3*h)^2/2; h^2/2 0 h^3/6 (2*h)^3/6 (3*h)^3/6; h^3/6 0 h^4/24 (2*h)^4/24 (3*h)^4/24];
B3 = [0;-1;0;0;0];

X3 = linsolve(A3,B3)

%% problem 4

kh = linspace(0,2*pi);
k_prime_square = 2-2*cos(kh);


figure
plot(kh,kh)
hold on
plot(kh,k_prime_square.^.5) 

% part b

A4 = [0 0 1 1 1; 0 0 -h 0 h; 1 1 h^2/2 0 h^2/2; -h h -h^3/6 0 h^3/6; h^2/2 h^2/2 h^4/24 0 h^4/24];