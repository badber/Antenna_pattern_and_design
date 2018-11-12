%% Short dipole

clear all
close all
clc

Il = 2;
r = 2;
lambda = 2e-3;

theta = 0:0.01:2*pi;
E_theta = abs(60*pi*Il/(r*lambda)*sin(theta));

figure(1)
polar(theta-pi/2,E_theta)% I substract pi/2 from theta so that the pattern corresponds to a vertical antenna instead of a horizontal one
title('E_{\theta} as a function of \theta')
