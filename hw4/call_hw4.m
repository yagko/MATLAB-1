clc 
clear all
close all

% Predefined test matrix (non-singular and well-defined) to test the function
A = [ 2  0 1
      3  2 5
      1 -1 0 ];

% Calculate the inverse of A using gaussian_inverter function
inv_A = gaussian_inverter( A )