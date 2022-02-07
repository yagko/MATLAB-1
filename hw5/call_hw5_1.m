clc;clear all;close all

% Specify the interval bounds
tmin = 5;
tmax = 15;

% Call the function
[ poly_deg , E_max ] = lagrange_int( tmin , tmax )