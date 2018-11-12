clear all
close all
clc

prompt = sprintf ('What is the number of elements in the array? \n');
n = input(prompt);
if ~isscalar(n) || (n~=floor(n)) || n<=1
    error('Introduce an integer scalar bigger than 1 for the number of elements in the array')
end

prompt = sprintf ('Broadside or end-fire array? \n(answer "broadside" or "end_fire", as a string in MATLAB) \n');
array_type = input(prompt);
if  ~strcmp(array_type,'broadside') && ~strcmp(array_type,'end_fire')
    error('Wrong type of array introduced.')
end

prompt = sprintf ('What is the value of "d", separation between \neach antenna and the following one in the array? \n(insert a number in terms of lambda, \nEx: introducing 1/2 means d=lambda/2) \n');
% If given kd or alpha, do this (see notes of lecture 18, page 4):
% If I am given the value of kd, calculate d as d=kd/2*pi.
% If I am given the value of alpha and it is different than 0, it means 
% that it is an end_fire array, then alpha=kd and d = alpha/2*pi.
d = input(prompt);
if ~isscalar(d)
    error('Introduce a scalar for the value of "d"')
end

kd = 2*pi*d;

if  strcmp(array_type,'broadside')
    alpha = 0; % Because alpha = 0 for broadside (Lecture 18 notes, page 1)
elseif strcmp(array_type,'end_fire')
    alpha = -kd; % Because alpha = -kd for end-fire (Lecture 18 notes, page 1)
end

%%
% Now calculate the radiation pattern as a function of phi:
n = (n+1)/2; % Because the "n" that I use in the expression for the triangular array is not actually its length, it is the length of the uniform array that it comes from (see my class notes back of page 20, the length of the triangular array is 2*n-1, being "n" the length of a uniform array that the triangular array is derived from) 
phi = 0:0.001:pi;
%psi = 0:0.0001:pi;
%psi = kd*cos(phi)+alpha;
%z = exp(j*psi);
z = exp(j*(kd*cos(phi)+alpha));
E = 0; % Initialize
for array_position=1:n
    E = E + array_position*z.^(array_position-1);
end
substract = 0; % Initialize. This is used to get the decreasing part of the triangle coefficients
for array_position=(n+1):(2*n-1)
    E = E + (n-1-substract)*z.^(array_position-1);
    substract = substract+1;
end

E = abs(E);

figure(1)
plot(phi,E)
%axis([0 pi 0 n])
xlabel('\phi','fontsize',15)
ylabel('E_T/E_o','fontsize',12)
title(['n=' num2str(2*n-1)])

% Now let's draw the polar plot:
% Now get the rest of the values for phi, since the radiation pattern as a
% function of phi is symmetrical with respect to the phi=0 axis:
l = length(phi);
phi(l+1:2*l) = -fliplr(phi);
E(l+1:2*l) = fliplr(E);

figure(2)
polar(phi,E)
title('E as a function of \phi')

clear i psi ET_over_E0

%%
% Calculate width of the primary maxima (Units: radians):
if  strcmp(array_type,'broadside')
    width_PM = 2/(n*d) % I don't include "lambda" in this formula because "d" has been specified in terms of "times of lambda", so it cancels out lambda
    disp('Units: radians')
    disp('NOTE: this used the formula for uniform distribution, I guess it is also valid for triangular')
elseif strcmp(array_type,'end_fire')
    width_PM = 2*sqrt(2/(n*d)) % I don't include "lambda" in this formula because "d" has been specified in terms of "times of lambda", so it cancels out lambda
    disp('Units: radians')
    disp('NOTE: this used the formula for uniform distribution, I guess it is also valid for triangular')
end
