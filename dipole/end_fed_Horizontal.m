% Author: Luis Badesa

%% Horizontal Dipole antenna (end-fed)
clear all
close all
clc

prompt = sprintf ('What is the length of the dipole antenna, in terms of lambda? \nL = x means here that L = x*lambda \n');
L = input(prompt); % "L = x" means here that "L = x*lambda"
if ~isscalar(L) || L<0
    error('Introduce a scalar bigger than 0 for the length of the dipole antenna')
end

% Modify these next 2 parameters for each case:
Im = 2;
r = 2;

% Dipole antenna equations:
theta = 0:0.01:2*pi;
u = 2*pi*L/2*(cos(theta)-1);
E_theta = abs(30*Im*sin(theta)*L.*sin(u)./u); % This formula is taken from page 5 of "traveling_wave.pdf"

figure(1)
polar(theta,E_theta)

% 3D plot:

[theta,phi] = meshgrid(0:0.01:pi,0:0.01:2*pi);
u = 2*pi*L/2*(cos(theta)-1);
E_theta = abs(30*Im*sin(theta)*L.*sin(u)./u);
 
% Transforming to XYZ (cartersian)
x = E_theta.*sin(theta).*cos(phi);
y = E_theta.*sin(theta).*sin(phi);
z = E_theta.*cos(theta);

% For cutting with 2 semi-planes:
a = x<0;
b = y<0;
c = a.*b;

for j=1:size(c,1)
    for k=1:size(c,2)
        if c(j,k)>0
            c(j,k) = NaN;
        end
    end
end
            
z_noPlot = z.*c;
z = z + z_noPlot;

figure(2)
mesh(x,y,z);
title('THIS IS FOR A VERTICAL ANTENNA')
