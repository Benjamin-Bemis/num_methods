%{

@author: Benjamin Bemis Ph.D Student, 
Advisor: Dr Juliano


Description:
AME 60614: Numerical Methods
Homework: 6
Due: 11/22/2024


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

%% Setting data paths
% Make sure to update this for the machine that you are working on. (Maybe, This should now run on any machine without change. 7/24/24)
% Change the current folder to the folder of this m-file.
if(~isdeployed)
  cd(fileparts(matlab.desktop.editor.getActiveFilename));
end

addpath(cd)
% cd ..; % Moving up a directory (from processing_code)
basepath = cd; % Pulling the current directory


imagepath = [basepath filesp 'images' filesep]; 
mkdir(imagepath);

   

%% Problem 21 Chapter 4
k1 = 0.04;
k2 = 10;
k3 = 1.5e3;

c1_0 = 0.9;
c2_0 = 0.1;
c3_0 = 0;



dc1 = @(c1,c2,c3) -k1*c1 +k2*c2*c3;
dc2 = @(c1,c2,c3) k1*c1 -k2*c2*c3 - 2*k3*c2^2;
dc3 = @(c1,c2,c3) 2*k3*c2^2;

% Sum of c1,c2,c3 is equal to 1






%% Problem 26 Chapter 4



%% Problem 27 Chapter 4